"""
F-actin / phalloidin readouts that do not depend on where the cell boundaries are drawn.

F-actin often sits between cells (stress fibres, junctions, aggregates, protrusions) where no
segmentation assigns it to a cell, and a mean intensity dilutes a few very bright structures. So for
every analysed channel whose stain is F-actin or phalloidin (with "Overall signal" quantified):

Per field (fields.csv; SuperPlots with one dot per field):
  <ch>_field_mean_bgsub                  mean background-subtracted intensity over the whole field
  <ch>_field_total_per_nucleus           total background-subtracted signal of the field / nuclei
  <ch>_field_bright_area_pct             % of the field above the bright cut-off
  <ch>_field_bright_area_per_nucleus_um2 bright area of the field / nuclei (µm²)
Per cell (nuclei.csv; fields with cell boundaries):
  <ch>_cell_bright_pct                   % of the cell's area above the bright cut-off
  <ch>_cell_p95_bgsub                    95th percentile of the cell's background-subtracted intensity

Bright cut-off: the `actin_bright_percentile` (default 99th) of the background-subtracted pixels inside
the reference condition's cells (e.g. WT · MOCK) of the same staining panel - "brighter than 99 % of the
reference cells' F-actin" - or a fixed value (`actin_bright_threshold`). Saturated pixels count as bright,
so the readouts hold even where the camera clipped.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import tifffile

ACTIN_RE = re.compile(r"actin(?!in)|phalloidin", re.IGNORECASE)  # F-actin, actin, phalloidin - not actinin
HIST_BINS = 4096


def is_actin(stain) -> bool:
    return bool(ACTIN_RE.search(str(stain or "")))


def actin_channels(fields, channels: list[str], params) -> list[str]:
    """Analysed channels stained for F-actin / phalloidin in at least one field, with the overall signal on."""
    out = []
    for c in channels:
        if "signal" in params.quantify_for(c) and any(is_actin(f.stains.get(c, "")) for f in fields):
            out.append(c)
    return out


def _edges(ceiling: float) -> np.ndarray:
    return np.linspace(0.0, float(ceiling), HIST_BINS + 1)


def reference_histogram(fi, nuclei_path: Path, cells_path: "Path | None", channel: str, px_um: float,
                        ceiling: float) -> np.ndarray:
    """Histogram (HIST_BINS bins over 0..ceiling) of one field's background-subtracted pixels of `channel`
    inside cells (or, without cell boundaries, over the whole field) - pooled over the reference fields."""
    from .measure import field_background

    nuclei = tifffile.imread(nuclei_path).astype(np.int32)
    cells = (tifffile.imread(cells_path).astype(np.int32)
             if cells_path is not None and Path(cells_path).exists() else None)
    img = tifffile.imread(fi.channels[channel])
    bg = field_background({channel: img}, nuclei, cells, px_um)[channel]
    v = img.astype(np.float64) - bg
    if cells is not None and (cells > 0).any():
        v = v[cells > 0]
    return np.histogram(np.clip(v, 0, ceiling), bins=_edges(ceiling))[0].astype(np.int64)


def threshold_from_histogram(hist: np.ndarray, ceiling: float, percentile: float) -> float:
    """The `percentile` of a HIST_BINS histogram (upper edge of the bin where it is reached)."""
    total = hist.sum()
    if total == 0:
        return float("nan")
    k = int(np.searchsorted(np.cumsum(hist), total * percentile / 100.0))
    return float(_edges(ceiling)[min(k + 1, HIST_BINS)])


def field_readouts(c: str, img: np.ndarray, bg: float, threshold: float, n_nuclei: int, px2: float) -> dict:
    """Whole-field F-actin readouts (no cell boundaries needed)."""
    v = np.clip(img.astype(np.float64) - bg, 0, None)
    bright = int((v > threshold).sum()) if np.isfinite(threshold) else 0
    nn = n_nuclei if n_nuclei > 0 else float("nan")
    return {f"{c}_field_mean_bgsub": float(v.mean()),
            f"{c}_field_total_per_nucleus": float(v.sum() / nn),
            f"{c}_field_bright_area_pct": 100.0 * bright / v.size if np.isfinite(threshold) else float("nan"),
            f"{c}_field_bright_area_per_nucleus_um2": bright * px2 / nn if np.isfinite(threshold) else float("nan"),
            f"{c}_bright_threshold": float(threshold)}


def cell_readouts(c: str, values_bgsub: np.ndarray, threshold: float) -> dict:
    """Per-cell bright-structure readouts from the cell's background-subtracted pixel values."""
    if not values_bgsub.size:
        return {f"{c}_cell_bright_pct": float("nan"), f"{c}_cell_p95_bgsub": float("nan")}
    return {f"{c}_cell_bright_pct": (100.0 * float((values_bgsub > threshold).mean())
                                     if np.isfinite(threshold) else float("nan")),
            f"{c}_cell_p95_bgsub": float(np.percentile(values_bgsub, 95))}
