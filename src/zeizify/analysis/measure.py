"""
Per-field measurement: given nuclear labels, detect spots in each signal channel and measure
every nucleus, every spot and every channel pair (colocalization).

Per nucleus and channel (prefix "<channel>_"):
  mean, median, integrated     raw intensity inside the nucleus
  bg, mean_bgsub               extranuclear background (10th percentile of pixels away from nuclei,
                               i.e. camera offset + autofluorescence) and the background-subtracted
                               mean - background-paired intensities as in SpeckleSPEC
  cv                           coefficient of variation of the background-subtracted nuclear signal -
                               a threshold-free dispersal readout: high when signal is concentrated
                               in bodies, low when it is diffuse
  enrichment                   mean (bg-subtracted) intensity in spots / in the rest of the nucleus
  frac_in_spots                fraction of the nuclear (bg-subtracted) signal that lies in spots
  spots_n, spots_area_um2, spots_area_frac, spot_mean_area_um2, spot_mean_diam_um,
  spot_mean_area_halfmax_um2, spot_mean_diam_halfmax_um  (spot size at half of its own peak - independent of
                               brightness, unlike the thresholded area; use it to compare sizes),
  spot_mean_ecc, spot_mean_solidity, spot_mean_intensity
                               spot counts, sizes and shapes (Alexander 2024; Aktas lab 2026)

Nucleus QC: nuclei whose mean DAPI is far from the field's median (condensed mitotic/apoptotic
or dim out-of-focus nuclei) or that contain saturated pixels are flagged (qc_pass = False). They
stay in the tables and are excluded from the statistics.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations
from pathlib import Path

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage import measure, morphology

from .coloc import costes_pvalue, object_coloc, pixel_coloc
from .io import FieldImages
from .params import AnalysisParams
from .segment import detect_spots, halfmax_areas


@dataclass
class FieldResult:
    key: str
    nuclei: list[dict] = field(default_factory=list)
    spots: list[dict] = field(default_factory=list)
    coloc: list[dict] = field(default_factory=list)
    field_row: dict = field(default_factory=dict)
    error: str = ""


def _ids(fi: FieldImages) -> dict:
    return {"well": fi.well, "field": fi.field, "antibody": fi.antibody, "treatment": fi.treatment,
            "cell_type": fi.cell_type, "replicate": fi.replicate, "label": fi.label}


def _nan_mean(x) -> float:
    x = np.asarray(x, dtype=float)
    return float(np.nanmean(x)) if x.size and np.isfinite(x).any() else float("nan")


def measure_field(fi: FieldImages, nuclei: np.ndarray, nuclear: str, spot_channels: list[str],
                  params: AnalysisParams, px_um: float, bit_ceiling: "int | None" = None,
                  qc_png: "Path | None" = None, spot_floors: "dict[str, float] | None" = None) -> FieldResult:
    res = FieldResult(fi.key)
    ids = _ids(fi)
    imgs = {c: tifffile.imread(fi.channels[c]) for c in [nuclear, *spot_channels] if c in fi.channels}
    spot_channels = [c for c in spot_channels if c in imgs]
    dapi = imgs[nuclear].astype(np.float64)
    def _ceiling(a: np.ndarray) -> float:
        if bit_ceiling:
            return bit_ceiling
        return np.iinfo(a.dtype).max if np.issubdtype(a.dtype, np.integer) else np.inf

    ceiling = {c: _ceiling(a) for c, a in imgs.items()}
    px2 = px_um * px_um

    # background from pixels well away from any nucleus
    away = ~morphology.dilation(nuclei > 0, morphology.disk(max(1, int(1.0 / px_um)), decomposition="sequence"))
    bg = {c: (float(np.percentile(a[away], 10)) if away.any() else float(np.percentile(a, 1))) for c, a in imgs.items()}

    spots, tophats = {}, {}
    for c in spot_channels:
        spots[c], tophats[c] = detect_spots(imgs[c], nuclei, params, px_um, (spot_floors or {}).get(c, 0.0))

    # spot table + per-nucleus spot summaries
    spot_summary: dict[tuple[str, int], list[dict]] = {}
    for c in spot_channels:
        lab = spots[c]
        if lab.max() == 0:
            continue
        parent = ndi.maximum(nuclei, lab, index=np.arange(1, lab.max() + 1))
        th = tophats[c]
        hm = halfmax_areas(lab, th, nuclei, int(params.spot_max_diameter_um / px_um / 2) + 2)
        for sp in measure.regionprops(lab, intensity_image=imgs[c]):
            nid = int(parent[sp.label - 1])
            coords = sp.coords
            row = {**ids, "channel": c, "stain": fi.stains.get(c, ""), "spot_id": sp.label, "nucleus_id": nid,
                   "area_um2": sp.area * px2, "eq_diam_um": sp.equivalent_diameter_area * px_um,
                   "area_halfmax_um2": hm.get(sp.label, sp.area) * px2,
                   "diam_halfmax_um": 2 * math.sqrt(hm.get(sp.label, sp.area) / math.pi) * px_um,
                   "eccentricity": sp.eccentricity, "solidity": sp.solidity,
                   "mean_intensity": float(sp.intensity_mean), "max_intensity": float(sp.intensity_max),
                   "integrated_intensity": float(sp.intensity_mean * sp.area),
                   "integrated_tophat": float(th[coords[:, 0], coords[:, 1]].sum()),
                   "centroid_y_um": sp.centroid[0] * px_um, "centroid_x_um": sp.centroid[1] * px_um}
            res.spots.append(row)
            spot_summary.setdefault((c, nid), []).append(row)

    # nuclei
    regions = measure.regionprops(nuclei)
    dapi_means = np.array([dapi[r.coords[:, 0], r.coords[:, 1]].mean() for r in regions]) if regions else np.array([])
    dapi_med = float(np.median(dapi_means)) if len(dapi_means) else float("nan")
    rng = np.random.default_rng(params.seed)
    block = max(1, int(round(params.psf_um / px_um)))
    pairs = list(combinations(spot_channels, 2))
    slices = ndi.find_objects(nuclei)

    for r, dmean in zip(regions, dapi_means):
        nid = r.label
        sl = slices[nid - 1]
        mask = nuclei[sl] == nid
        n = {**ids, "nucleus_id": nid, "centroid_y_um": r.centroid[0] * px_um, "centroid_x_um": r.centroid[1] * px_um,
             "area_um2": r.area * px2, "eccentricity": r.eccentricity, "solidity": r.solidity,
             "dapi_mean": float(dmean), "dapi_integrated": float(dmean * r.area),
             "dapi_ratio": float(dmean / dapi_med) if dapi_med else float("nan")}
        sat = max(float((imgs[c][sl][mask] >= ceiling[c]).mean()) for c in imgs)
        n["saturated_frac"] = sat
        reasons = []
        if n["dapi_ratio"] < params.qc_dapi_low:
            reasons.append("dim DAPI")
        if n["dapi_ratio"] > params.qc_dapi_high:
            reasons.append("bright DAPI (mitotic/apoptotic?)")
        if sat > params.qc_max_saturated_frac:
            reasons.append("saturated")
        n["qc_pass"] = not reasons
        n["qc_reason"] = "; ".join(reasons)

        for c in spot_channels:
            v = imgs[c][sl][mask].astype(np.float64)
            vb = v - bg[c]
            in_spot = spots[c][sl][mask] > 0
            rows = spot_summary.get((c, nid), [])
            mean_b = float(vb.mean())
            n.update({
                f"{c}_stain": fi.stains.get(c, ""),
                f"{c}_mean": float(v.mean()), f"{c}_median": float(np.median(v)), f"{c}_integrated": float(v.sum()),
                f"{c}_bg": bg[c], f"{c}_mean_bgsub": mean_b,
                f"{c}_cv": float(vb.std() / mean_b) if mean_b > 0 else float("nan"),
                f"{c}_spots_n": len(rows),
                f"{c}_spots_area_um2": float(in_spot.sum() * px2),
                f"{c}_spots_area_frac": float(in_spot.mean()),
                f"{c}_spot_mean_area_um2": _nan_mean([x["area_um2"] for x in rows]),
                f"{c}_spot_mean_diam_um": _nan_mean([x["eq_diam_um"] for x in rows]),
                f"{c}_spot_mean_area_halfmax_um2": _nan_mean([x["area_halfmax_um2"] for x in rows]),
                f"{c}_spot_mean_diam_halfmax_um": _nan_mean([x["diam_halfmax_um"] for x in rows]),
                f"{c}_spot_mean_ecc": _nan_mean([x["eccentricity"] for x in rows]),
                f"{c}_spot_mean_solidity": _nan_mean([x["solidity"] for x in rows]),
                f"{c}_spot_mean_intensity": _nan_mean([x["mean_intensity"] for x in rows]),
            })
            if in_spot.any() and (~in_spot).any():
                nucleo = float(vb[~in_spot].mean())
                n[f"{c}_enrichment"] = float(vb[in_spot].mean() / nucleo) if nucleo > 0 else float("nan")
            else:
                n[f"{c}_enrichment"] = float("nan")
            tot = float(np.clip(vb, 0, None).sum())
            n[f"{c}_frac_in_spots"] = float(np.clip(vb[in_spot], 0, None).sum() / tot) if tot > 0 else float("nan")
        res.nuclei.append(n)

        for a, b in pairs:
            # background-subtracted, clipped at 0, so Manders' "> 0" means "above background"
            va = np.clip(imgs[a][sl][mask].astype(np.float64) - bg[a], 0, None)
            vbb = np.clip(imgs[b][sl][mask].astype(np.float64) - bg[b], 0, None)
            row = {**ids, "nucleus_id": nid, "qc_pass": n["qc_pass"], "channel_a": a, "channel_b": b,
                   "stain_a": fi.stains.get(a, ""), "stain_b": fi.stains.get(b, ""),
                   "pair_type": _pair_type(fi.stains.get(a, ""), fi.stains.get(b, ""))}
            row.update(pixel_coloc(va, vbb))
            row["costes_p"] = (costes_pvalue(imgs[a][sl] - bg[a], imgs[b][sl] - bg[b], mask, block,
                                             params.coloc_randomizations, rng)
                               if params.coloc_randomizations > 0 else float("nan"))
            row.update(object_coloc(np.where(mask, spots[a][sl], 0), np.where(mask, spots[b][sl], 0), mask, px_um))
            res.coloc.append(row)

    res.field_row = {**ids, "nuclei_total": len(regions), "nuclei_qc_pass": sum(x["qc_pass"] for x in res.nuclei),
                     **{f"{c}_bg": bg[c] for c in imgs},
                     **{f"{c}_saturated_field": fi.saturated.get(c, False) for c in imgs}}

    if qc_png is not None:
        write_overlay(qc_png, imgs, nuclear, spot_channels, nuclei, spots, res.nuclei)
    return res


def _pair_type(sa: str, sb: str) -> str:
    if sa and sb:
        return "double-stained"
    if sa or sb:
        return "single-stain control"  # the unstained channel reports bleed-through / crosstalk
    return "unstained"


OVERLAY_COLORS = [(0, 255, 0), (255, 0, 255), (0, 255, 255), (255, 255, 0)]


def write_overlay(path: Path, imgs: dict, nuclear: str, spot_channels: list[str], nuclei: np.ndarray,
                  spots: dict, nuc_rows: list[dict], max_side: int = 1600) -> None:
    """Composite (nuclear stain blue, first two spot channels green/magenta) with nucleus outlines
    (white = analysed, red = QC-flagged) and spot outlines in each channel's color."""
    from PIL import Image
    from skimage import segmentation, transform

    def stretch(a):
        lo, hi = np.percentile(a, [0.5, 99.8])
        return np.clip((a.astype(np.float32) - lo) / max(hi - lo, 1), 0, 1)

    h, w = nuclei.shape
    s = min(1.0, max_side / max(h, w))
    shape = (max(1, int(h * s)), max(1, int(w * s)))

    def small(a, order=1):
        return transform.resize(a, shape, order=order, preserve_range=True, anti_aliasing=order > 0)

    rgb = np.zeros((*shape, 3), np.float32)
    rgb[..., 2] += 0.8 * small(stretch(imgs[nuclear]))
    for c, col in zip(spot_channels, OVERLAY_COLORS):
        rgb += 0.8 * small(stretch(imgs[c]))[..., None] * (np.array(col, np.float32) / 255)
    rgb = np.clip(rgb, 0, 1) * 255
    nl = small(nuclei, order=0).astype(np.int32)
    failed = {r["nucleus_id"] for r in nuc_rows if not r["qc_pass"]}
    edge = segmentation.find_boundaries(nl, mode="outer")
    near_failed = np.isin(ndi.grey_dilation(nl, size=3), list(failed)) if failed else np.zeros_like(edge)
    rgb[edge] = (255, 255, 255)
    rgb[edge & near_failed] = (255, 60, 60)
    for c, col in zip(spot_channels, OVERLAY_COLORS):
        se = segmentation.find_boundaries(small(spots[c], order=0).astype(np.int32), mode="inner")
        rgb[se] = col
    path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgb.astype(np.uint8)).save(path)
