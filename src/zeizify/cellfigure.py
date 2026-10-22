"""
Single-cell figures for publication: one nucleus per cell type x treatment, cropped from the
extracted TIFFs and laid out with the stain panels as rows or columns.

Choosing the cells
    Run Zeizify Analyze first: the measurements (analysis/nuclei.csv, colocalization.csv) and
    the nuclear masks (analysis/masks/) are what make a cell "representative". For each
    condition the candidates are ranked by how closely they match the condition's medians:
      * hard match: spot counts per channel and double-positive spots (+/- max(1, 15 %)), Pearson's r
        (+/- 0.1) and spot size at half maximum per channel (+/- 30 %); widened in steps only if
        nothing fits, and an isolated cell is preferred over a slightly closer match;
      * score: weighted robust distance (|value - median| / IQR) over Pearson's r, spot size,
        intensity CV and mean intensity per channel, nuclear area and DAPI ratio;
      * interphase: mean DAPI 0.75-1.25 x the field median (not condensed, not dim);
      * geometry: the whole nucleus inside the crop, few other nuclei in it, and at least one
        signal channel >= 3 x brighter inside the nucleus than outside (rejects rounded mitotic
        cells whose signal fills the cell).
    Any candidate can replace the automatic pick (manual choice).

Display
    Background-subtracted (field background: 10th percentile of pixels outside nuclei), linear,
    one range per channel shared by every condition (default) so panels are comparable, or one
    range per panel. Merges add the channels' pseudocolors.

Outputs (<out>/figures/<name>/): the figure (PNG / TIFF / PDF / SVG; editable text in PDF and
SVG), cells.csv (each cell and its values vs the condition medians), legend.txt (draft figure
legend) and figure_settings.json.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

from . import __version__
from .layout import COLORS, Panel, default_color_name, nuclear_channel

PAGE_WIDTHS_IN = {
    "Journal - 1 column (89 mm)": 89 / 25.4,
    "Journal - 1.5 columns (120 mm)": 120 / 25.4,
    "Journal - 2 columns (183 mm)": 183 / 25.4,
    "Slide 16:9 (13.33 in)": 13.333,
}
FORMATS = ("png", "tif", "pdf", "svg")
# (match tolerance x, max. fraction of the crop covered by other nuclei): an isolated cell is preferred over a
# slightly closer match
TIERS = ((1.0, 0.03), (1.5, 0.03), (1.0, 0.15), (2.0, 0.03), (1.5, 0.15), (2.0, 0.15), (3.0, 0.15))


# ---------------------------------------------------------------- data

@dataclass
class CellChoice:
    cell_type: str
    treatment: str
    well: str
    field: str
    nucleus_id: int
    how: str = "auto"  # "auto" or "manual"


@dataclass
class FigureOptions:
    panels: list[Panel]
    colors: dict[str, str] = field(default_factory=dict)
    cell_types: list[str] = field(default_factory=list)
    treatments: list[str] = field(default_factory=list)
    cell_type_labels: dict[str, str] = field(default_factory=dict)
    treatment_labels: dict[str, str] = field(default_factory=dict)
    crop_w_um: float = 30.0
    crop_h_um: float = 30.0
    rows: str = "treatment"  # "treatment": treatments as rows, panels as columns; "panel": the transpose
    blocks: str = "side"  # cell types side by side ("side") or stacked ("stacked")
    width_in: float = 183 / 25.4
    dpi: int = 300
    font: str = "Arial"
    font_size: float = 8.0
    gap_mm: float = 0.5
    block_gap_mm: float = 3.0
    contrast: str = "shared"  # "shared" per channel across conditions, or "panel"
    saturated_pct: float = 0.3  # % of nuclear pixels allowed to saturate when setting the maximum
    max_override: dict[str, float] = field(default_factory=dict)
    outline: str = "nuclear"  # outline the measured nucleus: "none", "nuclear" (nuclear-stain panels), "all"
    outline_color: str = "#FFD400"
    scale_bar_um: float = 5.0
    scale_bar: str = "last"  # "none", "last" (bottom-right panel) or "all"
    scale_bar_label: bool = True
    formats: list[str] = field(default_factory=lambda: ["png", "pdf"])

    def to_dict(self) -> dict:
        d = asdict(self)
        d["panels"] = ["+".join(p) for p in self.panels]
        return d


@dataclass
class CellData:
    out_dir: Path
    ds: object  # analysis.io.Dataset
    cells: pd.DataFrame  # QC-passed nuclei with colocalization columns
    signal: list[str]  # analysed signal channels
    nuclear: str

    def conditions(self) -> list[tuple[str, str]]:
        from .analysis.stats import is_control, is_wildtype, level_order

        ct = level_order(self.cells, "cell_type", first=is_wildtype)
        tr = level_order(self.cells, "treatment", first=is_control)
        have = set(zip(self.cells["cell_type"], self.cells["treatment"]))
        return [(c, t) for c in ct for t in tr if (c, t) in have]


def load_cells(out_dir: "str | Path") -> CellData:
    """Measurements of an analysed Zeizify output folder (Analyze must have been run)."""
    from .analysis.io import load_dataset

    out_dir = Path(out_dir)
    adir = out_dir / "analysis"
    if not (adir / "nuclei.csv").exists():
        raise FileNotFoundError("No analysis found - run Analyze on this output folder first "
                                f"({adir / 'nuclei.csv'} is missing).")
    ds = load_dataset(out_dir)
    if not ds.px_um:
        raise ValueError("Pixel size unknown - cannot crop in µm.")
    nuc = pd.read_csv(adir / "nuclei.csv")
    nuc = nuc[nuc["qc_pass"]].copy()
    for c in ("cell_type", "treatment", "well", "field"):
        nuc[c] = nuc[c].fillna("").astype(str)
    coloc_path = adir / "colocalization.csv"
    if coloc_path.exists():
        co = pd.read_csv(coloc_path)
        if len(co):
            first = co[["channel_a", "channel_b"]].drop_duplicates().iloc[0]
            co = co[(co["channel_a"] == first["channel_a"]) & (co["channel_b"] == first["channel_b"])]
            keep = [c for c in ("n_double", "pearson", "n_a_only", "n_b_only") if c in co]
            co = co[["well", "field", "nucleus_id", *keep]].astype({"well": str, "field": str})
            nuc = nuc.merge(co, on=["well", "field", "nucleus_id"], how="left")
    nuclear = ds.nuclear or ""
    signal = [c for c in ds.signal_channels if f"{c}_spots_n" in nuc]
    return CellData(out_dir, ds, nuc, signal, nuclear)


# ---------------------------------------------------------------- choosing cells

def _metrics(cd: CellData) -> tuple[dict[str, float], dict[str, float], list[str]]:
    """(hard count metrics -> tolerance kind, soft metrics -> weight, size metrics)."""
    cols = cd.cells.columns
    key = {f"{c}_spots_n": 1.0 for c in cd.signal}
    if "n_double" in cols:
        key["n_double"] = 1.0
    sizes = [f"{c}_spot_mean_area_halfmax_um2" for c in cd.signal if f"{c}_spot_mean_area_halfmax_um2" in cols]
    soft = {"pearson": 3.0, **{s: 2.0 for s in sizes}, **{f"{c}_cv": 1.0 for c in cd.signal},
            **{f"{c}_mean_bgsub": 1.0 for c in cd.signal}, "area_um2": 0.5, "dapi_ratio": 0.5}
    soft = {k: w for k, w in soft.items() if k in cols}
    return key, soft, sizes


class _FieldCache:
    """Masks and channel images of recently used fields (they are large)."""

    def __init__(self, cd: CellData, size: int = 4):
        self.cd, self.size, self._d = cd, size, {}

    def get(self, well: str, fld: str):
        k = (well, fld)
        if k in self._d:
            self._d[k] = self._d.pop(k)
            return self._d[k]
        fi = next(f for f in self.cd.ds.fields if f.well == well and f.field == fld)
        mask = tifffile.imread(self.cd.out_dir / "analysis" / "masks" / f"{fi.key}_nuclei.tif").astype(np.int32)
        imgs, bg = {}, {}
        outside = mask[::4, ::4] == 0
        for c in self.cd.ds.channels:
            if c in fi.channels:
                im = tifffile.imread(fi.channels[c]).astype(np.float32)
                imgs[c] = im
                sub = im[::4, ::4][outside]
                bg[c] = float(np.percentile(sub, 10)) if sub.size else float(np.percentile(im, 1))
        self._d[k] = (fi, mask, imgs, bg)
        while len(self._d) > self.size:
            self._d.pop(next(iter(self._d)))
        return self._d[k]


def _window(cy: int, cx: int, h: int, w: int) -> tuple[int, int]:
    return cy - h // 2, cx - w // 2


def _crop(arr: np.ndarray, y0: int, x0: int, h: int, w: int, fill=0) -> np.ndarray:
    out = np.full((h, w), fill, dtype=arr.dtype)
    ys, xs = max(0, y0), max(0, x0)
    ye, xe = min(arr.shape[0], y0 + h), min(arr.shape[1], x0 + w)
    if ye > ys and xe > xs:
        out[ys - y0:ye - y0, xs - x0:xe - x0] = arr[ys:ye, xs:xe]
    return out


def _geometry_ok(cd: CellData, cache: _FieldCache, r, h: int, w: int, neighbours_max: float) -> "tuple[bool, float]":
    fi, mask, imgs, bg = cache.get(r["well"], r["field"])
    px = cd.ds.px_um
    cy, cx = int(round(r["centroid_y_um"] / px)), int(round(r["centroid_x_um"] / px))
    y0, x0 = _window(cy, cx, h, w)
    if y0 < 0 or x0 < 0 or y0 + h > mask.shape[0] or x0 + w > mask.shape[1]:
        return False, 1.0
    crop = mask[y0:y0 + h, x0:x0 + w]
    own = crop == int(r["nucleus_id"])
    if not own.any():
        return False, 1.0
    other = float(((crop > 0) & ~own).mean())
    m = 0.04 * min(h, w)
    rr, cc = np.nonzero(own)
    if not (rr.min() > m and cc.min() > m and rr.max() < h - m and cc.max() < w - m) or other >= neighbours_max:
        return False, other
    enriched = False
    for c in cd.signal:
        sig = np.clip(imgs[c][y0:y0 + h, x0:x0 + w] - bg[c], 0, None)
        outside = sig[crop == 0]
        if outside.size and sig[own].mean() >= 3 * max(float(outside.mean()), 1e-9):
            enriched = True
            break
    return enriched, other


def rank_candidates(cd: CellData, cell_type: str, treatment: str, crop_w_um: float, crop_h_um: float,
                    n: int = 12, cache: "_FieldCache | None" = None) -> pd.DataFrame:
    """Up to n candidate nuclei for one condition, best first, with their key values and the
    condition medians (columns '<metric>' and 'median_<metric>'), the tolerance tier that admitted
    them ('relax') and the fraction of the crop covered by other nuclei ('neighbours')."""
    cache = cache or _FieldCache(cd)
    cond = cd.cells[(cd.cells["cell_type"] == cell_type) & (cd.cells["treatment"] == treatment)]
    key, soft, sizes = _metrics(cd)
    meds = {c: float(np.nanmedian(cond[c])) if cond[c].notna().any() else float("nan") for c in [*key, *soft]}
    iqr = {c: (float(np.nanpercentile(cond[c], 75) - np.nanpercentile(cond[c], 25)) if cond[c].notna().any() else 1.0)
           or 1.0 for c in soft}
    px = cd.ds.px_um
    h, w = int(round(crop_h_um / px)), int(round(crop_w_um / px))
    base = cond[cond["dapi_ratio"].between(0.75, 1.25)] if "dapi_ratio" in cond else cond
    score = np.zeros(len(base))
    for c, wt in soft.items():
        score += wt * np.nan_to_num(np.abs(base[c].to_numpy(dtype=float) - meds[c]) / iqr[c], nan=3.0)
    base = base.assign(score=score)
    picked, seen = [], set()
    for relax, neigh in TIERS:
        d = base
        for c in key:
            if np.isfinite(meds[c]):
                d = d[(d[c] - meds[c]).abs() <= relax * max(1.0, 0.15 * meds[c])]
        for s in sizes:
            if np.isfinite(meds[s]):
                d = d[d[s].isna() | ((d[s] - meds[s]).abs() <= relax * 0.3 * meds[s])]
        if "pearson" in meds and np.isfinite(meds["pearson"]):
            d = d[(d["pearson"] - meds["pearson"]).abs() <= relax * 0.1]
        for _, r in d.sort_values("score").head(200).iterrows():
            ident = (r["well"], r["field"], int(r["nucleus_id"]))
            if ident in seen:
                continue
            ok, other = _geometry_ok(cd, cache, r, h, w, neigh)
            if ok:
                seen.add(ident)
                picked.append({**r.to_dict(), "relax": relax, "neighbours": round(other, 3)})
                if len(picked) >= n:
                    break
        if len(picked) >= n:
            break
    out = pd.DataFrame(picked)
    for c, v in meds.items():
        out[f"median_{c}"] = v
    return out


def describe(row, cd: CellData) -> str:
    """One-line summary of a candidate vs its condition medians."""
    parts = [f"{row['well']} {row['field']} #{int(row['nucleus_id'])}"]
    for c in cd.signal:
        k = f"{c}_spots_n"
        if k in row and np.isfinite(row.get(f"median_{k}", np.nan)):
            parts.append(f"{c} {int(row[k])} ({row[f'median_{k}']:.0f})")
    if "n_double" in row and np.isfinite(row.get("median_n_double", np.nan)):
        parts.append(f"double {int(row['n_double'])} ({row['median_n_double']:.0f})")
    if "pearson" in row and np.isfinite(row.get("pearson", np.nan)):
        parts.append(f"r {row['pearson']:.2f} ({row.get('median_pearson', np.nan):.2f})")
    return "  ·  ".join(parts)


def auto_choices(cd: CellData, conditions: list[tuple[str, str]], crop_w_um: float, crop_h_um: float,
                 n: int = 12, progress=None) -> dict[tuple[str, str], pd.DataFrame]:
    """Ranked candidates for every condition (the first row is the automatic pick)."""
    cache = _FieldCache(cd)
    out = {}
    for i, (ct, t) in enumerate(conditions, start=1):
        out[(ct, t)] = rank_candidates(cd, ct, t, crop_w_um, crop_h_um, n, cache)
        if progress:
            progress(i, len(conditions), f"{ct} · {t}")
    return out


def choice_from_row(row, how: str = "auto") -> CellChoice:
    return CellChoice(str(row["cell_type"]), str(row["treatment"]), str(row["well"]), str(row["field"]),
                      int(row["nucleus_id"]), how)


# ---------------------------------------------------------------- rendering

def _font(name: str) -> str:
    from matplotlib import font_manager

    have = {f.name for f in font_manager.fontManager.ttflist}
    for cand in (name, "Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"):
        if cand in have:
            return cand
    return "DejaVu Sans"


def _caption(spec: Panel, stains: dict[str, str], nuclear: str, short: bool = False) -> str:
    if len(spec) == 1:
        c = spec[0]
        return "DAPI" if c == nuclear and not stains.get(c) else (f"{stains[c]}\n({c})" if stains.get(c) else c)
    if short:  # column header: the channels are named in their own columns
        return "Merge"
    names = [stains.get(c) or c for c in spec]
    return "Merge\n(" + " + ".join(names) + ")"


def _fit(fig, text, width_in: float) -> None:
    """Shrink a label's font until it fits width_in inches (never below 5 pt)."""
    w = text.get_window_extent(renderer=fig.canvas.get_renderer()).width / fig.dpi
    if w > width_in:
        text.set_fontsize(max(5.0, text.get_fontsize() * width_in / w))


def _hex(rgb) -> str:
    return "#{:02x}{:02x}{:02x}".format(*(int(round(v * 255)) for v in rgb))


def _label_color(spec: Panel, colors: dict[str, str]) -> str:
    if len(spec) != 1:
        return "#000000"
    rgb = COLORS[colors.get(spec[0]) or default_color_name(spec[0])]
    if rgb == (1.0, 1.0, 1.0):
        return "#404040"
    # darker shade of the channel color so it reads on white
    return _hex(tuple(0.72 * v for v in rgb))


def crops_for(cd: CellData, choices: dict[tuple[str, str], CellChoice], opts: FigureOptions,
              cache: "_FieldCache | None" = None) -> dict:
    """{(ct, t): {"img": {channel: bg-subtracted float32}, "own": bool mask, "choice": CellChoice}}."""
    px = cd.ds.px_um
    h, w = int(round(opts.crop_h_um / px)), int(round(opts.crop_w_um / px))
    cache = cache or _FieldCache(cd)
    chans = sorted({c for p in opts.panels for c in p})
    out = {}
    for k, ch in choices.items():
        fi, mask, imgs, bg = cache.get(ch.well, ch.field)
        row = cd.cells[(cd.cells["well"] == ch.well) & (cd.cells["field"] == ch.field)
                       & (cd.cells["nucleus_id"] == ch.nucleus_id)]
        if row.empty:  # e.g. QC-failed: centroid from the mask
            yy, xx = np.nonzero(mask == ch.nucleus_id)
            cy, cx = int(yy.mean()), int(xx.mean())
        else:
            cy, cx = int(round(row["centroid_y_um"].iloc[0] / px)), int(round(row["centroid_x_um"].iloc[0] / px))
        y0, x0 = _window(cy, cx, h, w)
        out[k] = {"img": {c: np.clip(_crop(imgs[c], y0, x0, h, w) - bg[c], 0, None) for c in chans if c in imgs},
                  "own": _crop(mask, y0, x0, h, w) == ch.nucleus_id, "choice": ch}
    return out


def thumbnail(cd: CellData, choice: CellChoice, opts: FigureOptions, size: int = 96,
              cache: "_FieldCache | None" = None):
    """Small PIL preview of one cell: the merge of the figure's signal channels (or the first panel),
    scaled to its own nucleus, with the nucleus outlined."""
    from PIL import Image
    from scipy import ndimage as ndi

    merge = next((p for p in opts.panels if len(p) > 1), None) or (opts.panels[0] if opts.panels else ())
    one = FigureOptions(panels=[merge], colors=opts.colors, crop_w_um=opts.crop_w_um, crop_h_um=opts.crop_h_um,
                        contrast="panel", saturated_pct=opts.saturated_pct)
    k = (choice.cell_type, choice.treatment)
    crop = crops_for(cd, {k: choice}, one, cache)[k]
    # one range per channel for every condition (as in the figure), from the measured nuclear means: a dim,
    # diffuse condition stays dim instead of being stretched
    rng = {}
    for c in merge:
        col = f"{c}_mean_bgsub"
        if col in cd.cells and cd.cells[col].notna().any():
            rng[c] = max(4.0 * float(cd.cells.groupby(["cell_type", "treatment"])[col].median().max()), 1.0)
        else:
            rng.update({cc: v for (*_, cc), v in display_ranges({k: crop}, one).items() if cc == c})
    arr = (_rgb(crop, merge, rng, k, opts.colors) * 255).astype(np.uint8)
    own = crop["own"]
    arr[own & ~ndi.binary_erosion(own, iterations=3)] = (255, 212, 0)
    img = Image.fromarray(arr)
    img.thumbnail((size, size))
    return img


def display_ranges(crops: dict, opts: FigureOptions) -> dict:
    """{channel: max} shared across conditions, or {(ct, t, channel): max} per panel."""
    q = 100.0 - opts.saturated_pct
    chans = sorted({c for v in crops.values() for c in v["img"]})
    if opts.contrast == "panel":
        out = {}
        for k, v in crops.items():
            for c in chans:
                px = v["img"][c][v["own"]] if v["own"].any() else v["img"][c].ravel()
                out[(*k, c)] = opts.max_override.get(c) or max(float(np.percentile(px, q)), 1.0)
        return out
    out = {}
    for c in chans:
        pooled = np.concatenate([v["img"][c][v["own"]] if v["own"].any() else v["img"][c].ravel() for v in crops.values()])
        out[c] = opts.max_override.get(c) or max(float(np.percentile(pooled, q)), 1.0)
    return out


def _rgb(crop: dict, spec: Panel, rng: dict, key: tuple, colors: dict[str, str]) -> np.ndarray:
    acc = np.zeros((*crop["own"].shape, 3), dtype=np.float32)
    for c in spec:
        vmax = rng.get(c) or rng.get((*key, c)) or 1.0
        acc += np.clip(crop["img"][c] / vmax, 0, 1)[..., None] * np.array(COLORS[colors.get(c) or default_color_name(c)],
                                                                             dtype=np.float32)
    return np.clip(acc, 0, 1)


def render(cd: CellData, choices: dict[tuple[str, str], CellChoice], opts: FigureOptions,
           base: "Path | None" = None, preview_dpi: "int | None" = None):
    """Draw the figure. With base (a path without extension) the figure is saved in every format of
    opts.formats and the list of files is returned; otherwise the matplotlib Figure is returned."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.transforms import offset_copy

    family = _font(opts.font)
    plt.rcParams.update({"font.family": family, "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fs = opts.font_size
    cts = [c for c in opts.cell_types if any(k[0] == c for k in choices)]
    trs = [t for t in opts.treatments if any(k[1] == t for k in choices)]
    panels = [p for p in opts.panels if p]
    crops = crops_for(cd, {k: v for k, v in choices.items() if k[0] in cts and k[1] in trs}, opts)
    rng = display_ranges(crops, opts)
    stains = {}
    for f in cd.ds.fields:
        stains.update({c: s for c, s in f.stains.items() if s})
    tlabel = lambda t: opts.treatment_labels.get(t, t)  # noqa: E731
    clabel = lambda c: opts.cell_type_labels.get(c, c)  # noqa: E731

    # grid inside one block
    if opts.rows == "panel":
        R, C = len(panels), len(trs)
        row_lab = [(_caption(p, stains, cd.nuclear), _label_color(p, opts.colors)) for p in panels]
        col_lab = [(tlabel(t), "#000000") for t in trs]
    else:
        R, C = len(trs), len(panels)
        row_lab = [(tlabel(t), "#000000") for t in trs]
        col_lab = [(_caption(p, stains, cd.nuclear, short=True), _label_color(p, opts.colors)) for p in panels]
    B = max(1, len(cts))
    side = opts.blocks == "side" or B == 1
    pt = 1 / 72
    gap, bgap = opts.gap_mm / 25.4, opts.block_gap_mm / 25.4
    lines_col = max(t.count("\n") + 1 for t, _ in col_lab) if col_lab else 1
    head_h = lines_col * fs * 1.25 * pt + 4 * pt
    title_h = (fs * 1.3 * 1.5 * pt + 6 * pt) if (B > 1 or any(cts)) else 0.0
    longest = max((max(len(s) for s in t.split("\n")) for t, _ in row_lab), default=4)
    left = longest * fs * 0.56 * pt + 8 * pt
    right, bottom = 4 * pt, 4 * pt
    aspect = opts.crop_h_um / opts.crop_w_um
    ncols_total = B * C if side else C
    avail = opts.width_in - left - right - (ncols_total - B if side else C - 1) * gap - ((B - 1) * bgap if side else 0)
    p = avail / ncols_total
    if p <= 0.05:
        raise ValueError("Figure too narrow for this many panels - increase the width or the number of columns.")
    ph = p * aspect
    block_h = R * ph + (R - 1) * gap
    if side:
        H = title_h + head_h + block_h + bottom
    else:
        H = head_h + B * (title_h + block_h) + (B - 1) * bgap + bottom
    W = opts.width_in
    fig = plt.figure(figsize=(W, H), facecolor="white")

    def ax_at(x, y_top, w, h):
        return fig.add_axes([x / W, 1 - (y_top + h) / H, w / W, h / H])

    last_ax = None
    all_axes = []
    for bi, ct in enumerate(cts or [""]):
        if side:
            x0 = left + bi * (C * p + (C - 1) * gap + bgap)
            ty, y0 = 0.0, title_h + head_h
        else:  # stacked: [title][column headers - first block only][panels]
            x0 = left
            ty = 0.0 if bi == 0 else title_h + head_h + block_h + bgap + (bi - 1) * (title_h + block_h + bgap)
            y0 = ty + title_h + (head_h if bi == 0 else 0.0)
        bw = C * p + (C - 1) * gap
        if title_h:
            fig.text((x0 + bw / 2) / W, 1 - (ty + title_h * 0.42) / H, clabel(ct), ha="center", va="center",
                     fontsize=fs * 1.3, fontweight="bold")
            fig.add_artist(plt.Line2D([x0 / W, (x0 + bw) / W], [1 - (ty + title_h * 0.85) / H] * 2, color="black",
                                      lw=0.8))
        if bi == 0 or side:  # column headers (once on top when stacked)
            hy = title_h + head_h - 2 * pt
            for ci, (txt, col) in enumerate(col_lab):
                cxp = x0 + ci * (p + gap) + p / 2
                _fit(fig, fig.text(cxp / W, 1 - hy / H, txt, ha="center", va="bottom", fontsize=fs, color=col,
                                   fontweight="bold", linespacing=1.1), p * 0.98)
        for ri in range(R):
            for ci in range(C):
                t = trs[ri] if opts.rows != "panel" else trs[ci]
                spec = panels[ci] if opts.rows != "panel" else panels[ri]
                ax = ax_at(x0 + ci * (p + gap), y0 + ri * (ph + gap), p, ph)
                ax.set_xticks([])
                ax.set_yticks([])
                for s in ax.spines.values():
                    s.set_visible(False)
                k = (ct, t)
                if k not in crops or not all(c in crops[k]["img"] for c in spec):
                    ax.set_facecolor("#f2f2f2")
                    ax.text(0.5, 0.5, "n/a", ha="center", va="center", fontsize=fs, color="0.5", transform=ax.transAxes)
                    continue
                ax.imshow(_rgb(crops[k], spec, rng, k, opts.colors), interpolation="nearest")
                is_nuc = spec == (cd.nuclear,)
                if opts.outline == "all" or (opts.outline == "nuclear" and is_nuc):
                    ax.contour(crops[k]["own"].astype(float), levels=[0.5], colors=opts.outline_color,
                               linewidths=0.6, linestyles="--")
                all_axes.append(ax)
                last_ax = ax
                if ci == 0 and (bi == 0 or not side):
                    txt, col = row_lab[ri]
                    fig.text((x0 - 5 * pt) / W, 1 - (y0 + ri * (ph + gap) + ph / 2) / H, txt, ha="right", va="center",
                             fontsize=fs, color=col, fontweight="bold", linespacing=1.1)

    if opts.scale_bar != "none" and opts.scale_bar_um > 0 and all_axes:
        hpx, wpx = next(iter(crops.values()))["own"].shape
        n = opts.scale_bar_um / cd.ds.px_um
        for ax in (all_axes if opts.scale_bar == "all" else [last_ax]):
            xr = wpx - n - 0.06 * wpx
            yb = hpx * 0.93
            ax.plot([xr, xr + n], [yb, yb], color="white", lw=max(1.5, fs * 0.28), solid_capstyle="butt")
            if opts.scale_bar_label:
                tr = offset_copy(ax.transData, fig=fig, y=2, units="points")
                ax.text(xr + n / 2, yb, f"{opts.scale_bar_um:g} µm", color="white", ha="center", va="bottom",
                        fontsize=fs * 0.85, transform=tr)
            ax.set_xlim(-0.5, wpx - 0.5)
            ax.set_ylim(hpx - 0.5, -0.5)

    if base is None:
        if preview_dpi:
            fig.set_dpi(preview_dpi)
        return fig, crops, rng
    base = Path(base)
    base.parent.mkdir(parents=True, exist_ok=True)
    files = []
    for fmt in opts.formats:
        path = base.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=opts.dpi, facecolor="white")
        if fmt in ("png", "tif"):  # flatten to RGB (journals reject alpha) and record the resolution
            from PIL import Image

            with Image.open(path) as im:
                rgb = im.convert("RGB")
            extra = {"compression": "tiff_lzw"} if fmt == "tif" else {}
            rgb.save(path, dpi=(opts.dpi, opts.dpi), **extra)
        files.append(path)
    plt.close(fig)
    return files, crops, rng


# ---------------------------------------------------------------- export

def export(cd: CellData, choices: dict[tuple[str, str], CellChoice], opts: FigureOptions, name: str = "cell_figure",
           candidates: "dict | None" = None) -> dict:
    """Write the figure, cells.csv, legend.txt and figure_settings.json to <out>/figures/<name>/."""
    from .plate import slug

    folder = cd.out_dir / "figures" / (slug(name) or "cell_figure")
    files, crops, rng = render(cd, choices, opts, folder / (slug(name) or "cell_figure"))
    rows = []
    key, soft, sizes = _metrics(cd)
    for (ct, t), ch in choices.items():
        if (ct, t) not in crops:
            continue
        cond = cd.cells[(cd.cells["cell_type"] == ct) & (cd.cells["treatment"] == t)]
        row = cd.cells[(cd.cells["well"] == ch.well) & (cd.cells["field"] == ch.field)
                       & (cd.cells["nucleus_id"] == ch.nucleus_id)]
        rec = {"cell_type": ct, "treatment": t, "well": ch.well, "field": ch.field, "nucleus_id": ch.nucleus_id,
               "chosen": ch.how, "nuclei_in_condition": len(cond)}
        for c in [*key, *sizes, "pearson"]:
            if c in cd.cells:
                rec[c] = float(row[c].iloc[0]) if len(row) and pd.notna(row[c].iloc[0]) else float("nan")
                rec[f"median_{c}"] = float(np.nanmedian(cond[c])) if cond[c].notna().any() else float("nan")
        rows.append(rec)
    pd.DataFrame(rows).to_csv(folder / "cells.csv", index=False)
    (folder / "legend.txt").write_text(legend_text(cd, opts, rng), encoding="utf-8")
    settings = {"created": datetime.now(timezone.utc).isoformat(timespec="seconds"), "zeizify": __version__,
                "input": str(cd.out_dir), "pixel_size_um": cd.ds.px_um, "options": opts.to_dict(),
                "cells": [asdict(c) for c in choices.values()],
                "display_max": {("|".join(k) if isinstance(k, tuple) else k): v for k, v in rng.items()}}
    (folder / "figure_settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")
    return {"folder": folder, "files": files}


def legend_text(cd: CellData, opts: FigureOptions, rng: dict) -> str:
    stains = {}
    for f in cd.ds.fields:
        stains.update({c: s for c, s in f.stains.items() if s})
    col = lambda c: (opts.colors.get(c) or default_color_name(c)).lower()  # noqa: E731
    chans = sorted({c for p in opts.panels for c in p}, key=lambda c: (c == cd.nuclear, c))
    names = [f"{stains.get(c) or ('DNA' if c == cd.nuclear else c)} ({c}, {col(c)})" for c in chans]
    shared = opts.contrast == "shared"
    txt = (
        "Representative single cells. "
        + ("Each panel shows one nucleus per cell type and treatment, chosen automatically from the QC-passed, "
           "interphase nuclei of that condition as the one closest to the condition's median spot counts, "
           "double-positive spots, spot size and Pearson's r (Zeizify Analyze). ")
        + f"Channels: {', '.join(names)}. Merges are additive overlays of the pseudocolored channels. "
        + "Images are maximum-intensity projections, background-subtracted and linearly scaled; "
        + ("the display range of each channel is identical in all conditions. " if shared
           else "each panel is scaled individually, so intensities are not comparable between panels. ")
        + ("Dashed outline: the nucleus measured. " if opts.outline != "none" else "")
        + (f"Scale bar, {opts.scale_bar_um:g} µm. " if opts.scale_bar != "none" else "")
        + f"Crops are {opts.crop_w_um:g} × {opts.crop_h_um:g} µm."
    )
    ranges = "; ".join(f"{k}: 0-{v:.0f}" for k, v in rng.items() if isinstance(k, str))
    return txt + ("\n\nDisplay ranges (background-subtracted grey levels): " + ranges if ranges else "") + "\n"


def default_figure_panels(channels: list[str]) -> list[Panel]:
    """Nuclear stain, each signal channel, then all signal channels merged."""
    nuc = nuclear_channel(channels)
    sig = [c for c in channels if c != nuc]
    out: list[Panel] = ([(nuc,)] if nuc else []) + [(c,) for c in sig]
    if len(sig) >= 2:
        out.append(tuple(sig))
    return out


def parse_order(text: str) -> list[Panel]:
    """'DAPI | AF488 | AF488+AF594' -> [('DAPI',), ('AF488',), ('AF488', 'AF594')]."""
    from .layout import parse_panel

    return [parse_panel(t) for t in text.replace(",", "|").split("|") if t.strip()]


def parse_labels(text: str) -> dict[str, str]:
    """'ML=ML792; MG=MG132' -> {'ML': 'ML792', 'MG': 'MG132'}."""
    out = {}
    for part in text.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            if k.strip():
                out[k.strip()] = v.strip()
    return out


def mm(x_in: float) -> float:
    return x_in * 25.4


__all__ = ["CellChoice", "FigureOptions", "PAGE_WIDTHS_IN", "auto_choices", "choice_from_row", "default_figure_panels",
           "describe", "export", "load_cells", "parse_labels", "parse_order", "rank_candidates", "render",
           "thumbnail"]
