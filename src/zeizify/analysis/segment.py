"""
Segmentation: nuclei from DAPI, then spots (speckles / nuclear bodies) inside each nucleus.

Nuclei
  cellpose   Cellpose generalist model (Stringer et al., Nat Methods 2021; 'nuclei' model in v2/v3,
             'cpsam' in v4), run on an image downscaled so nuclei are ~30 px across.
  classical  The CellProfiler IdentifyPrimaryObjects recipe: Gaussian smoothing, Otsu threshold,
             hole filling, then declumping of touching nuclei by a distance-transform watershed.
             Boundaries are refined at full resolution.
  Both are followed by a size window (fraction of the typical nuclear area) and removal of nuclei
  touching the image edge (Aktas lab speckle workflow, Cell 2026).

Spots - the "nuclear mask, then speckle objects" logic of the speckle papers:
  1. per nucleus, white top-hat background subtraction (structuring disk larger than the largest
     spot), which removes the smooth nucleoplasmic signal;
  2. threshold = max(Otsu of the top-hat values in that nucleus, median + k * 1.4826 * MAD).
     Per-object thresholding adapts to cell-to-cell staining differences; the MAD floor keeps
     dispersed nuclei (no real speckles, Ilik et al. eLife 2020) from being split into noise spots;
  3. optional declumping of touching spots at intensity saddles (watershed);
  4. size window on equivalent diameter, in µm.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation, transform
from skimage.feature import peak_local_max

from .params import AnalysisParams

_CELLPOSE_MODEL = None


def cellpose_available() -> bool:
    try:
        import cellpose  # noqa: F401

        return True
    except Exception:  # noqa: BLE001
        return False


def resolve_method(params: AnalysisParams) -> str:
    if params.nucleus_method == "auto":
        return "cellpose" if cellpose_available() else "classical"
    return params.nucleus_method


def _scale_for(diameter_px: float, target_px: float) -> float:
    return min(1.0, target_px / max(diameter_px, 1.0))


def _upscale_labels(lab: np.ndarray, shape: tuple[int, int]) -> np.ndarray:
    if lab.shape == shape:
        return lab
    return transform.resize(lab, shape, order=0, preserve_range=True, anti_aliasing=False).astype(np.int32)


def _filter_nuclei(lab: np.ndarray, params: AnalysisParams, px_um: float) -> np.ndarray:
    if params.exclude_border_nuclei:
        lab = segmentation.clear_border(lab)
    expected = params.nucleus_area_um2() / (px_um ** 2)
    lo, hi = params.nucleus_min_area_frac * expected, params.nucleus_max_area_frac * expected
    areas = np.bincount(lab.ravel())
    bad = np.where((areas < lo) | (areas > hi))[0]
    bad = bad[bad > 0]
    if len(bad):
        lab = lab.copy()
        lab[np.isin(lab, bad)] = 0
    lab, _, _ = segmentation.relabel_sequential(lab)
    return lab.astype(np.int32)


def segment_nuclei_classical(dapi: np.ndarray, params: AnalysisParams, px_um: float) -> np.ndarray:
    d_px = params.nucleus_diameter_um / px_um
    s = _scale_for(d_px, 40.0)
    img = dapi.astype(np.float32)
    small = transform.rescale(img, s, anti_aliasing=True, preserve_range=True) if s < 1 else img
    sm = filters.gaussian(small, sigma=1.0, preserve_range=True)
    t = filters.threshold_otsu(sm)
    mask = ndi.binary_fill_holes(sm > t)
    d_small = d_px * s
    mask = morphology.remove_small_objects(mask, max_size=max(4, int(0.05 * math.pi * (d_small / 2) ** 2)))
    dist = ndi.distance_transform_edt(mask)
    dist = filters.gaussian(dist, sigma=max(1.0, d_small / 20), preserve_range=True)
    cc = measure.label(mask)
    peaks = peak_local_max(dist, min_distance=max(2, int(0.3 * d_small)), labels=cc, exclude_border=False)
    markers = np.zeros(mask.shape, dtype=np.int32)
    markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
    lab = segmentation.watershed(-dist, markers, mask=mask)
    lab = _upscale_labels(lab, dapi.shape)
    if s < 1:  # refine the blocky upscaled boundary with the full-resolution foreground
        full = filters.gaussian(img, sigma=1.0 / s * 0.5, preserve_range=True) > t
        lab[~full] = 0
    return _filter_nuclei(lab, params, px_um)


def segment_nuclei_cellpose(dapi: np.ndarray, params: AnalysisParams, px_um: float) -> np.ndarray:
    global _CELLPOSE_MODEL
    import cellpose
    from cellpose import models

    major = int(str(getattr(cellpose, "__version__", None) or _cellpose_version()).split(".")[0])
    d_px = params.nucleus_diameter_um / px_um
    s = _scale_for(d_px, 30.0)
    img = dapi.astype(np.float32)
    small = transform.rescale(img, s, anti_aliasing=True, preserve_range=True) if s < 1 else img
    if _CELLPOSE_MODEL is None:
        _CELLPOSE_MODEL = (models.CellposeModel(gpu=params.use_gpu) if major >= 4
                           else models.CellposeModel(gpu=params.use_gpu, model_type="nuclei"))
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # torch sparse-tensor notices from Cellpose's dynamics
        if major >= 4:
            out = _CELLPOSE_MODEL.eval(small, diameter=d_px * s)
        else:
            out = _CELLPOSE_MODEL.eval(small, diameter=d_px * s, channels=[0, 0])
    lab = np.asarray(out[0]).astype(np.int32)
    return _filter_nuclei(_upscale_labels(lab, dapi.shape), params, px_um)


def _cellpose_version() -> str:
    from importlib.metadata import version

    return version("cellpose")


def segment_nuclei(dapi: np.ndarray, params: AnalysisParams, px_um: float, method: "str | None" = None) -> np.ndarray:
    method = method or resolve_method(params)
    if method == "cellpose":
        return segment_nuclei_cellpose(dapi, params, px_um)
    return segment_nuclei_classical(dapi, params, px_um)


def detect_spots(img: np.ndarray, nuclei: np.ndarray, params: AnalysisParams, px_um: float,
                 min_brightness: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Spot labels (unique across the field, 0 = background) and the white-top-hat image used to
    find them (background-subtracted signal, zero outside nuclei). `min_brightness` > 0 also drops spots whose
    mean top-hat intensity is below it (the fixed, experiment-wide threshold)."""
    min_d = params.spot_min_diameter_um / px_um
    max_d = params.spot_max_diameter_um / px_um
    radius = max(2, int(round(0.75 * max_d)))  # disk wider than the largest spot's radius
    fp = morphology.disk(radius, decomposition="sequence")
    sigma = max(0.5, min_d / 4)
    spots = np.zeros(img.shape, dtype=np.int32)
    tophat = np.zeros(img.shape, dtype=np.float32)
    next_id = 1
    min_area, max_area = math.pi * (min_d / 2) ** 2, math.pi * (max_d / 2) ** 2

    for reg in measure.regionprops(nuclei):
        r0, c0, r1, c1 = reg.bbox
        R0, C0 = max(0, r0 - radius), max(0, c0 - radius)
        R1, C1 = min(img.shape[0], r1 + radius), min(img.shape[1], c1 + radius)
        crop = img[R0:R1, C0:C1].astype(np.float32)
        th = morphology.white_tophat(filters.gaussian(crop, sigma=sigma, preserve_range=True), fp)
        inside = nuclei[R0:R1, C0:C1] == reg.label
        vals = th[inside]
        if vals.size < 10:
            continue
        med = float(np.median(vals))
        mad = 1.4826 * float(np.median(np.abs(vals - med)))
        floor = med + params.spot_sensitivity * mad
        try:
            otsu = float(filters.threshold_otsu(vals))
        except ValueError:
            otsu = floor
        thr = max(otsu, floor)
        fg = (th > thr) & inside
        if not fg.any():
            continue
        if params.spot_declump:
            pk = peak_local_max(np.where(fg, th, 0), min_distance=max(1, int(round(min_d / 2))),
                                labels=measure.label(fg), exclude_border=False)
            mk = np.zeros(fg.shape, dtype=np.int32)
            mk[tuple(pk.T)] = np.arange(1, len(pk) + 1)
            lab = segmentation.watershed(-th, mk, mask=fg)
        else:
            lab = measure.label(fg)
        areas = np.bincount(lab.ravel())
        keep = np.where((areas >= min_area) & (areas <= max_area))[0]
        keep = keep[keep > 0]
        if min_brightness > 0 and len(keep):
            sums = np.bincount(lab.ravel(), weights=th.ravel(), minlength=len(areas))
            keep = keep[sums[keep] / areas[keep] >= min_brightness]
        if len(keep):
            lut = np.zeros(lab.max() + 1, dtype=np.int32)
            lut[keep] = np.arange(next_id, next_id + len(keep), dtype=np.int32)
            new = lut[lab]
            view = spots[R0:R1, C0:C1]
            view[new > 0] = new[new > 0]
            next_id += len(keep)
        tophat[R0:R1, C0:C1][inside] = th[inside]
    return spots, tophat


def spot_brightness(spots: np.ndarray, tophat: np.ndarray) -> np.ndarray:
    """Mean top-hat intensity of every spot in a label image."""
    if spots.max() == 0:
        return np.zeros(0)
    area = np.bincount(spots.ravel())
    total = np.bincount(spots.ravel(), weights=tophat.ravel())
    ids = np.nonzero(area[1:])[0] + 1
    return total[ids] / area[ids]


SEPARABILITY_MIN = 0.7


def auto_spot_threshold(brightness: np.ndarray) -> tuple[float, float]:
    """Experiment-wide spot brightness floor: Otsu's threshold on log10 spot brightness pooled over all
    wells, and Otsu's separability (between-class / total variance, 0-1). A single bell-shaped
    population scores ~0.64 (2/pi for a normal distribution) and two distinct populations (dim grains vs
    bodies) score higher; below SEPARABILITY_MIN a fixed floor should not be used."""
    v = np.log10(brightness[brightness > 0])
    if v.size < 10 or np.ptp(v) == 0:
        return 0.0, float("nan")
    t = float(filters.threshold_otsu(v))
    lo, hi = v[v <= t], v[v > t]
    w0, w1 = lo.size / v.size, hi.size / v.size
    eta = float(w0 * w1 * (lo.mean() - hi.mean()) ** 2 / v.var()) if lo.size and hi.size else float("nan")
    return float(10 ** t), eta
