"""
Colocalization of two channels inside one nucleus.

Pixel-based - definitions follow Fiji's Coloc 2 (github.com/fiji/Colocalisation_Analysis):
  pearson          Pearson's r over all nuclear pixels.
  spearman         Spearman's rank correlation.
  icq              Li's intensity correlation quotient: fraction of pixels with
                   (a - mean a)(b - mean b) > 0, minus 0.5 (range -0.5 .. 0.5).
  costes thresholds  orthogonal (total-least-squares) regression b = m a + c; the threshold on the
                   stepped channel (a if |m| <= 1, else b; the other follows the line) is lowered by
                   bisection until Pearson's r of the pixels below threshold - a < tA OR b < tB,
                   strict - reaches 0 (Costes et al., Biophys J 2004).
  m1, m2           Manders: sum of a where b > 0 / sum of a (and vice versa). Zeizify passes
                   background-subtracted intensities clipped at 0, so "> 0" means above background;
                   inside a nucleus this is still close to 1 - prefer tm1 / tm2.
  tm1, tm2         thresholded Manders: sum of a where b >= tB (and b > 0) / sum of a; likewise tM2.
  costes_p         optional significance: the fraction of block-scrambled images (blocks ~ PSF size)
                   whose r is >= the observed r (Costes 2004); reported as 1 - P-value in Coloc 2.

Object-based - in the spirit of SpotitPy (BMC Bioinformatics 2022) and DiAna (Methods 2017):
  frac_a_overlap   fraction of channel-A spots overlapping at least one channel-B spot (and B->A).
  obj_m1, obj_m2   area of A spots inside B spots / area of A spots (and vice versa).
  enrich_a         obj_m1 divided by the fraction of the nucleus covered by B spots - the overlap
                   expected if A spots were placed at random; > 1 means more overlap than chance.
  nn_ab_um         median distance from each A spot centroid to the nearest B spot centroid.
"""

from __future__ import annotations

import math

import numpy as np
from scipy import ndimage as ndi
from scipy import stats
from scipy.spatial import cKDTree


def _pearson(a: np.ndarray, b: np.ndarray) -> float:
    if a.size < 3:
        return float("nan")
    da, db = a - a.mean(), b - b.mean()
    den = math.sqrt(float((da * da).sum()) * float((db * db).sum()))
    return float((da * db).sum() / den) if den > 0 else float("nan")


def regression_tls(a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    """Orthogonal regression b = m a + c, as in Coloc 2's AutoThresholdRegression."""
    va, vb = float(a.var()), float(b.var())
    cov = float(((a - a.mean()) * (b - b.mean())).mean())
    if cov == 0:
        return float("nan"), float("nan")
    m = (vb - va + math.sqrt((vb - va) ** 2 + 4 * cov * cov)) / (2 * cov)
    return m, float(b.mean() - m * a.mean())


def costes_thresholds(a: np.ndarray, b: np.ndarray, iters: int = 40) -> tuple[float, float]:
    """(tA, tB), or (nan, nan) when the channels are not positively correlated."""
    m, c = regression_tls(a, b)
    if not np.isfinite(m) or m <= 0 or _pearson(a, b) <= 0:
        return float("nan"), float("nan")
    step_a = abs(m) <= 1
    x = a if step_a else b
    lo, hi = float(x.min()), float(x.max())

    def thr(t):
        return (t, m * t + c) if step_a else ((t - c) / m, t)

    for _ in range(iters):
        mid = (lo + hi) / 2
        ta, tb = thr(mid)
        below = (a < ta) | (b < tb)
        r = _pearson(a[below], b[below]) if below.sum() >= 3 else float("nan")
        if not np.isfinite(r) or r <= 0:
            lo = mid
        else:
            hi = mid
    return thr((lo + hi) / 2)


def pixel_coloc(a: np.ndarray, b: np.ndarray) -> dict[str, float]:
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    out = {"pearson": _pearson(a, b)}
    if a.size >= 3 and a.std() > 0 and b.std() > 0:
        out["spearman"] = float(stats.spearmanr(a, b).statistic)
    else:
        out["spearman"] = float("nan")
    da, db = a - a.mean(), b - b.mean()
    out["icq"] = float(((da * db) > 0).mean() - 0.5) if a.size else float("nan")
    sa, sb = a.sum(), b.sum()
    out["m1"] = float(a[b > 0].sum() / sa) if sa > 0 else float("nan")
    out["m2"] = float(b[a > 0].sum() / sb) if sb > 0 else float("nan")
    ta, tb = costes_thresholds(a, b)
    out["costes_ta"], out["costes_tb"] = ta, tb
    if np.isfinite(ta):
        out["tm1"] = float(a[(b > 0) & (b >= tb)].sum() / sa) if sa > 0 else float("nan")
        out["tm2"] = float(b[(a > 0) & (a >= ta)].sum() / sb) if sb > 0 else float("nan")
    else:
        out["tm1"] = out["tm2"] = float("nan")
    return out


def costes_pvalue(a2d: np.ndarray, b2d: np.ndarray, mask: np.ndarray, block: int, n: int,
                  rng: np.random.Generator) -> float:
    """Fraction of block-scrambled A images with r >= observed (one-sided). Blocks are PSF-sized
    squares of the nucleus bounding box; the nuclear mask is applied after scrambling."""
    if n <= 0:
        return float("nan")
    a2d, b2d, mask = a2d.astype(np.float64), b2d.astype(np.float64), mask.astype(bool)
    h, w = a2d.shape
    block = max(1, int(block))
    H, W = (h // block) * block, (w // block) * block
    if H < block or W < block:
        return float("nan")
    a2d, b2d, mask = a2d[:H, :W], b2d[:H, :W], mask[:H, :W]
    obs = _pearson(a2d[mask], b2d[mask])
    blocks = a2d.reshape(H // block, block, W // block, block).swapaxes(1, 2).reshape(-1, block, block)
    ge = 0
    for _ in range(n):
        shuf = blocks[rng.permutation(len(blocks))]
        img = shuf.reshape(H // block, W // block, block, block).swapaxes(1, 2).reshape(H, W)
        ge += _pearson(img[mask], b2d[mask]) >= obs
    return float(ge / n)


def object_coloc(la: np.ndarray, lb: np.ndarray, nucleus: np.ndarray, px_um: float) -> dict[str, float]:
    """la, lb: spot label crops (0 = background); nucleus: boolean mask of the same crop."""
    ids_a = np.unique(la[la > 0])
    ids_b = np.unique(lb[lb > 0])
    out = {"n_a": len(ids_a), "n_b": len(ids_b)}
    in_a, in_b = la > 0, lb > 0
    both = in_a & in_b
    out["n_a_overlap"] = int(len(np.unique(la[both])))  # A spots touching any B spot
    out["n_b_overlap"] = int(len(np.unique(lb[both])))
    # double-positive spots: connected A-or-B regions containing both, so one PML body overlapping two
    # SUMO foci counts once
    union, _ = ndi.label(in_a | in_b)
    out["n_double"] = int(len(np.unique(union[both])))
    out["frac_a_overlap"] = float(out["n_a_overlap"] / len(ids_a)) if len(ids_a) else float("nan")
    out["frac_b_overlap"] = float(out["n_b_overlap"] / len(ids_b)) if len(ids_b) else float("nan")
    area_a, area_b, area_n = int(in_a.sum()), int(in_b.sum()), int(nucleus.sum())
    out["obj_m1"] = float(both.sum() / area_a) if area_a else float("nan")
    out["obj_m2"] = float(both.sum() / area_b) if area_b else float("nan")
    fb, fa = (area_b / area_n if area_n else 0), (area_a / area_n if area_n else 0)
    out["enrich_a"] = float(out["obj_m1"] / fb) if area_a and fb > 0 else float("nan")
    out["enrich_b"] = float(out["obj_m2"] / fa) if area_b and fa > 0 else float("nan")
    if len(ids_a) and len(ids_b):
        ca = np.array(_centroids(la, ids_a))
        cb = np.array(_centroids(lb, ids_b))
        out["nn_ab_um"] = float(np.median(cKDTree(cb).query(ca)[0]) * px_um)
        out["nn_ba_um"] = float(np.median(cKDTree(ca).query(cb)[0]) * px_um)
    else:
        out["nn_ab_um"] = out["nn_ba_um"] = float("nan")
    return out


def _centroids(lab: np.ndarray, ids: np.ndarray) -> list[tuple[float, float]]:
    from scipy import ndimage as ndi

    return ndi.center_of_mass(np.ones_like(lab, dtype=np.uint8), lab, ids)
