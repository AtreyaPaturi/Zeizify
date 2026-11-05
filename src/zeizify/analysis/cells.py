"""
Whole cells and cytoplasm from a cell-body / edge marker (F-actin, phalloidin, CellMask, tubulin, WGA)
plus the nuclear stain.

  cellpose   Cellpose-SAM ('cpsam', Cellpose 4; Pachitariu et al. 2025), a generalist model that
             separates touching cells, given [marker, nuclear stain]; 'cyto3' in Cellpose 2/3.
             The image is downscaled so a typical cell is ~90 px across.
  classical  seeded watershed from the nuclei on the marker's gradient, within the marker's
             foreground and one typical cell diameter of a nucleus (CellProfiler
             IdentifySecondaryObjects-like).

Each cell is tied to its nucleus: cell id = nucleus id. A cell that holds several nuclei (touching
cells found as one) is split between them by nearest nucleus; a cell without an analysed nucleus is
dropped. Cytoplasm = cell minus nuclei. Cells cut by the image edge keep only their nucleus, so no
cytoplasm is measured for them.
"""

from __future__ import annotations

import warnings

import numpy as np
from scipy import ndimage as ndi
from skimage import filters, segmentation, transform

from .params import AnalysisParams
from .segment import _cellpose_version, _scale_for, _upscale_labels, cellpose_available

# stain names that mark the cell body or edge; a channel whose antibody / dye matches gives the cell boundaries
CELL_MARKER_KEYWORDS = ("actin", "phalloidin", "cellmask", "cell mask", "tubulin", "wga", "membrane", "cytoplasm",
                        "cytosol")
_MODEL = None


def cell_marker_for(stains: dict[str, str], channels, params: AnalysisParams) -> "str | None":
    """The channel used for cell boundaries in a field: params.cell_marker if set, else the channel whose
    stain names a cell-body / edge marker. None = no cell segmentation for this field."""
    if params.cell_segmentation == "off":
        return None
    if params.cell_marker:
        return params.cell_marker if params.cell_marker in channels else None
    for c in channels:
        s = (stains.get(c) or "").lower()
        if s and any(k in s for k in CELL_MARKER_KEYWORDS):
            return c
    return None


def resolve_cell_method(params: AnalysisParams) -> str:
    if params.cell_method == "auto":
        return "cellpose" if cellpose_available() else "classical"
    return params.cell_method


def assign_cells(cells: np.ndarray, nuclei: np.ndarray) -> np.ndarray:
    """Relabel a whole-cell label image so each cell carries the id of its nucleus (see module notes)."""
    out = np.zeros(nuclei.shape, dtype=np.int32)
    if cells.max() > 0 and nuclei.max() > 0:
        base = int(nuclei.max()) + 1
        n_area = np.bincount(nuclei.ravel(), minlength=base)
        both = (cells > 0) & (nuclei > 0)
        ids, counts = np.unique(cells[both].astype(np.int64) * base + nuclei[both], return_counts=True)
        members: dict[int, list[int]] = {}
        for code, k in zip(ids, counts):
            c, nu = int(code // base), int(code % base)
            if k >= 0.5 * n_area[nu]:  # the nucleus lies (mostly) in this cell
                members.setdefault(c, []).append(nu)
        slices = ndi.find_objects(cells)
        for c, nus in members.items():
            sl = slices[c - 1]
            region = cells[sl] == c
            if len(nus) == 1:
                out[sl][region] = nus[0]
                continue
            seeds = np.where(np.isin(nuclei[sl], nus), nuclei[sl], 0)
            _, (iy, ix) = ndi.distance_transform_edt(seeds == 0, return_indices=True)
            out[sl][region] = seeds[iy, ix][region]
    out[nuclei > 0] = nuclei[nuclei > 0]
    return out


def drop_border_cells(cells: np.ndarray, nuclei: np.ndarray) -> np.ndarray:
    """Cells cut by the image edge keep only their nucleus."""
    edge = np.unique(np.concatenate([cells[0], cells[-1], cells[:, 0], cells[:, -1]]))
    edge = edge[edge > 0]
    if len(edge):
        cells = cells.copy()
        cells[np.isin(cells, edge) & (nuclei == 0)] = 0
    return cells


def segment_cells_cellpose(marker: np.ndarray, dapi: np.ndarray, params: AnalysisParams, px_um: float) -> np.ndarray:
    global _MODEL
    import cellpose
    from cellpose import models

    from . import segment

    major = int(str(getattr(cellpose, "__version__", None) or _cellpose_version()).split(".")[0])
    d_px = params.cell_diameter_um / px_um
    s = _scale_for(d_px, 90.0)
    stack = np.stack([transform.rescale(a.astype(np.float32), s, anti_aliasing=True, preserve_range=True) if s < 1
                      else a.astype(np.float32) for a in (marker, dapi)], -1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # torch sparse-tensor notices from Cellpose's dynamics
        if major >= 4:  # one generalist model (shared with nuclear segmentation)
            if segment._CELLPOSE_MODEL is None:
                segment._CELLPOSE_MODEL = models.CellposeModel(gpu=params.use_gpu)
            out = segment._CELLPOSE_MODEL.eval(stack, channel_axis=-1, diameter=None)
        else:
            if _MODEL is None:
                _MODEL = models.CellposeModel(gpu=params.use_gpu, model_type="cyto3")
            out = _MODEL.eval(stack, channel_axis=-1, channels=[1, 2], diameter=d_px * s)
    return _upscale_labels(np.asarray(out[0]).astype(np.int32), marker.shape)


def segment_cells_classical(marker: np.ndarray, nuclei: np.ndarray, params: AnalysisParams, px_um: float) -> np.ndarray:
    d_px = params.cell_diameter_um / px_um
    s = _scale_for(d_px, 90.0)
    img = marker.astype(np.float32)
    small = transform.rescale(img, s, anti_aliasing=True, preserve_range=True) if s < 1 else img
    nuc = transform.resize(nuclei, small.shape, order=0, preserve_range=True, anti_aliasing=False).astype(np.int32)
    sm = filters.gaussian(small, sigma=1.0, preserve_range=True)
    fg = ndi.binary_fill_holes(sm > filters.threshold_otsu(sm)) | (nuc > 0)
    fg &= ndi.distance_transform_edt(nuc == 0) <= d_px * s
    lab = segmentation.watershed(filters.sobel(sm), markers=nuc, mask=fg)
    return _upscale_labels(lab.astype(np.int32), marker.shape)


def segment_cells(marker: np.ndarray, dapi: np.ndarray, nuclei: np.ndarray, params: AnalysisParams, px_um: float,
                  method: "str | None" = None) -> np.ndarray:
    """Cell label image aligned to the nuclei (cell id = nucleus id)."""
    method = method or resolve_cell_method(params)
    raw = (segment_cells_cellpose(marker, dapi, params, px_um) if method == "cellpose"
           else segment_cells_classical(marker, nuclei, params, px_um))
    cells = assign_cells(raw, nuclei)
    return drop_border_cells(cells, nuclei) if params.exclude_border_cells else cells
