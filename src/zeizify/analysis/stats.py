"""
Statistics on replicate-level summaries, and SuperPlots.

Following Lord et al. (J Cell Biol 2020, "SuperPlots"), nuclei are not independent samples: every
metric is first summarised per replicate (median over that replicate's QC-passed nuclei), and the
tests compare those replicate medians. By default a replicate is a well - wells on one plate are
technical replicates; add a "Replicate" column to the plate map for biological replicates.

  omnibus    Welch's one-way ANOVA across treatments (unequal variances; scipy f_oneway).
  post hoc   Games-Howell pairwise comparisons (unequal variances and group sizes), as used with
             CellProfiler output by Kono et al. 2023 (Kimura lab).
  effect     Hedges' g on replicate medians, and Cliff's delta on all nuclei (descriptive only).
  multiple   Benjamini-Hochberg q-values across all omnibus tests of a run.
Treatments are compared within each staining panel and cell type, and cell types (e.g. WT vs KO)
within each staining panel and treatment, so only identically stained wells are compared.
"""

from __future__ import annotations

import math
import re
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import stats

from ..contrast import CONTROL_KEYWORDS

NUCLEUS_METRICS = {"area_um2": "nuclear area (µm²)", "dapi_integrated": "integrated DAPI"}
CHANNEL_METRICS = {
    "spots_n": "spots per nucleus",
    "spots_area_frac": "nuclear area in spots (fraction)",
    "spot_mean_area_um2": "mean spot area, thresholded (µm²)",
    "spot_mean_area_halfmax_um2": "mean spot size at half max (µm²)",
    "spot_mean_ecc": "mean spot eccentricity",
    "cv": "intensity CV in nucleus (dispersal)",
    "enrichment": "spot / nucleoplasm intensity",
    "frac_in_spots": "fraction of signal in spots",
    "mean_bgsub": "mean nuclear intensity (bg-subtracted)",
}
COLOC_METRICS = {
    "pearson": "Pearson's r",
    "spearman": "Spearman's ρ",
    "icq": "Li's ICQ",
    "tm1": "Manders tM1 (Costes)",
    "tm2": "Manders tM2 (Costes)",
    "n_double": "double-positive spots per nucleus",
    "n_a_overlap": "A spots overlapping B (count)",
    "n_b_overlap": "B spots overlapping A (count)",
    "frac_a_overlap": "A spots overlapping B (fraction)",
    "frac_b_overlap": "B spots overlapping A (fraction)",
    "enrich_a": "A-in-B overlap vs chance",
    "nn_ab_um": "A→B nearest spot (µm)",
}
SUPERPLOT_CHANNEL = ["spots_n", "spot_mean_area_halfmax_um2", "spots_area_frac", "cv", "enrichment"]
SUPERPLOT_COLOC = ["n_double", "pearson", "tm1", "tm2", "frac_a_overlap"]
WILDTYPE_KEYWORDS = ("wt", "wildtype", "wild", "parental", "control", "ctrl")
REPLICATE_COLORS = ["#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#0072B2", "#D55E00", "#F0E442", "#999999"]  # Okabe-Ito


def is_control(treatment) -> bool:
    return any(k in str(treatment).lower() for k in CONTROL_KEYWORDS)


def is_wildtype(cell_type) -> bool:
    return any(t in WILDTYPE_KEYWORDS for t in re.findall(r"[a-z0-9]+", str(cell_type).lower()))


def _well_key(well) -> tuple:
    m = re.fullmatch(r"([A-Za-z]+)0*(\d+)", str(well).strip())
    return (0, m.group(1).upper(), int(m.group(2))) if m else (1, str(well), 0)


def level_order(df: pd.DataFrame, col: str, preferred=(), first=None) -> list[str]:
    """Levels of `col` in display order: `preferred` (case-insensitive, in that order), then levels for
    which `first(level)` is true (controls / wild type), then plate-map order (the first well, row by
    row, in which the level appears) - e.g. MOCK, ML792, MG132, ML792+MG132 down a plate column."""
    if col not in df or not len(df):
        return []
    levels = list(dict.fromkeys(df[col].astype(str)))
    pos = {}
    if "well" in df:
        for lv, w in df.groupby(df[col].astype(str))["well"]:
            pos[lv] = min(_well_key(x) for x in w.astype(str).unique())
    pref = [str(p).strip().lower() for p in preferred if str(p).strip()]

    def key(lv):
        lo = lv.lower()
        return (pref.index(lo) if lo in pref else len(pref), 0 if first and first(lv) else 1, pos.get(lv, (2, lv, 0)))

    return sorted(levels, key=key)


def default_orders(long: pd.DataFrame, treatments=(), cell_types=()) -> dict:
    """Display / comparison order of treatments (control first) and cell types (wild type first)."""
    return {"treatment": level_order(long, "treatment", treatments, is_control),
            "cell_type": level_order(long, "cell_type", cell_types, is_wildtype)}


def long_table(nuclei: pd.DataFrame, coloc: pd.DataFrame, channels: list[str]) -> pd.DataFrame:
    """One row per nucleus x metric: antibody, cell_type, treatment, replicate, well, measure,
    metric, feature (channel or pair), stain, value. QC-failed nuclei are dropped."""
    ids = ["antibody", "cell_type", "treatment", "replicate", "well"]
    parts = []
    n = nuclei[nuclei["qc_pass"]] if len(nuclei) else nuclei
    for m, label in NUCLEUS_METRICS.items():
        if m in n:
            parts.append(n[ids].assign(metric=m, metric_label=label, feature="nucleus", stain="", value=n[m]))
    for c in channels:
        for m, label in CHANNEL_METRICS.items():
            col = f"{c}_{m}"
            if col in n:
                stain = n.get(f"{c}_stain", pd.Series("", index=n.index)).fillna("")
                parts.append(n[ids].assign(metric=m, metric_label=label, feature=c, stain=stain.values, value=n[col]))
    if len(coloc):
        cl = coloc[coloc["qc_pass"]]
        for m, label in COLOC_METRICS.items():
            if m in cl:
                feat = cl["channel_a"] + "×" + cl["channel_b"]
                stain = cl["stain_a"].fillna("").astype(str) + "×" + cl["stain_b"].fillna("").astype(str)
                parts.append(cl[ids].assign(metric=m, metric_label=label, feature=feat.values, stain=stain.values,
                                            value=cl[m]))
    if not parts:
        return pd.DataFrame(columns=ids + ["metric", "metric_label", "feature", "stain", "value"])
    out = pd.concat(parts, ignore_index=True)
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    return out


GROUP = ["antibody", "cell_type", "metric", "metric_label", "feature"]
OMNI_COLUMNS = GROUP + ["stain", "n_treatments", "replicates_per_treatment", "welch_F", "welch_p", "note"]
PAIR_COLUMNS = GROUP + ["treatment", "compare", "stain", "group_a", "group_b", "vs_control", "replicates_a",
                        "replicates_b", "nuclei_a", "nuclei_b", "median_nuclei_a", "median_nuclei_b",
                        "diff_replicate_means", "hedges_g_replicates", "cliffs_delta_nuclei", "gh_t", "gh_df", "gh_p",
                        "note"]


def replicate_summary(long: pd.DataFrame) -> pd.DataFrame:
    g = long.groupby(GROUP + ["stain", "treatment", "replicate"], dropna=False)["value"]
    rep = g.agg(n_nuclei="count", median="median", mean="mean").reset_index()
    return rep[rep["n_nuclei"] > 0]


def condition_summary(rep: pd.DataFrame) -> pd.DataFrame:
    g = rep.groupby(GROUP + ["stain", "treatment"], dropna=False)
    out = g.agg(n_replicates=("median", "count"), mean_of_medians=("median", "mean"),
                sd_of_medians=("median", "std"), n_nuclei=("n_nuclei", "sum")).reset_index()
    return out


def welch_anova(groups: list[np.ndarray]) -> tuple[float, float]:
    res = stats.f_oneway(*groups, equal_var=False)
    return float(res.statistic), float(res.pvalue)


def games_howell(groups: dict[str, np.ndarray]) -> list[dict]:
    """Pairwise Games-Howell; p from the studentized range distribution with k = number of groups."""
    k = len(groups)
    out = []
    for (na, a), (nb, b) in combinations(groups.items(), 2):
        ma, mb = a.mean(), b.mean()
        va, vb = a.var(ddof=1), b.var(ddof=1)
        sa, sb = va / len(a), vb / len(b)
        se = math.sqrt(sa + sb)
        if se == 0 or not np.isfinite(se):
            out.append({"group_a": na, "group_b": nb, "diff": mb - ma, "se": float("nan"), "t": float("nan"),
                        "df": float("nan"), "p": float("nan")})
            continue
        t = (mb - ma) / se
        df = (sa + sb) ** 2 / (sa ** 2 / (len(a) - 1) + sb ** 2 / (len(b) - 1))
        p = float(stats.studentized_range.sf(abs(t) * math.sqrt(2), k, df))
        out.append({"group_a": na, "group_b": nb, "diff": mb - ma, "se": se, "t": t, "df": df, "p": min(1.0, p)})
    return out


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float("nan")
    sp = math.sqrt(((na - 1) * a.var(ddof=1) + (nb - 1) * b.var(ddof=1)) / (na + nb - 2))
    if sp == 0:
        return float("nan")
    j = 1 - 3 / (4 * (na + nb) - 9)
    return float((b.mean() - a.mean()) / sp * j)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """P(b > a) - P(b < a), computed from ranks (O(n log n))."""
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if not len(a) or not len(b):
        return float("nan")
    u = stats.mannwhitneyu(b, a, alternative="two-sided").statistic
    return float(2 * u / (len(a) * len(b)) - 1)


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if ok.any():
        pv = p[ok]
        order = np.argsort(pv)
        ranked = pv[order] * len(pv) / np.arange(1, len(pv) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        out = np.empty_like(pv)
        out[order] = np.clip(ranked, 0, 1)
        q[ok] = out
    return q


def _pair_rows(base: dict, groups: dict, sel: pd.DataFrame, col: str, ref: str) -> list[dict]:
    """Pairwise rows for `groups` (level of `col` -> replicate medians): descriptive effect sizes always,
    Games-Howell only when every group has >= 2 replicates."""
    testable = {t: v for t, v in groups.items() if len(v) >= 2}
    gh = ({(r["group_a"], r["group_b"]): r for r in games_howell(testable)}
          if len(testable) >= 2 and len(testable) == len(groups) else {})
    rows = []
    for a, b in combinations(list(groups), 2):
        va = sel.loc[sel[col] == a, "value"].to_numpy(dtype=float)
        vb = sel.loc[sel[col] == b, "value"].to_numpy(dtype=float)
        g = gh.get((a, b), {})
        rows.append({**base, "group_a": a, "group_b": b, "vs_control": bool(ref) and ref in (a, b),
                     "replicates_a": len(groups[a]), "replicates_b": len(groups[b]),
                     "nuclei_a": int(np.isfinite(va).sum()), "nuclei_b": int(np.isfinite(vb).sum()),
                     "median_nuclei_a": float(np.nanmedian(va)) if np.isfinite(va).any() else float("nan"),
                     "median_nuclei_b": float(np.nanmedian(vb)) if np.isfinite(vb).any() else float("nan"),
                     "diff_replicate_means": float(groups[b].mean() - groups[a].mean()),
                     "hedges_g_replicates": hedges_g(groups[a], groups[b]),
                     "cliffs_delta_nuclei": cliffs_delta(va, vb),
                     "gh_t": g.get("t", float("nan")), "gh_df": g.get("df", float("nan")),
                     "gh_p": g.get("p", float("nan")),
                     "note": "" if g else "descriptive only: needs >= 2 replicates per group for a test"})
    return rows


def _select(df: pd.DataFrame, cols, vals) -> pd.DataFrame:
    for c, v in zip(cols, vals):
        df = df[df[c] == v]
    return df


def _ordered(levels: set, order: list) -> list:
    return [x for x in order if x in levels] + sorted(levels - set(order))


def run_tests(rep: pd.DataFrame, long: pd.DataFrame, orders: "dict | None" = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(omnibus table, pairwise table). Pairwise rows compare treatments within a cell type
    (compare = "treatment") and cell types within a treatment (compare = "cell_type")."""
    orders = orders or default_orders(long)
    omni, pairs = [], []
    for key, sub in rep.groupby(GROUP, dropna=False):
        order = _ordered(set(sub["treatment"]), orders["treatment"])
        groups = {t: sub.loc[sub["treatment"] == t, "median"].dropna().to_numpy() for t in order}
        groups = {t: v for t, v in groups.items() if len(v)}
        base = dict(zip(GROUP, key))
        base["stain"] = ";".join(sorted(set(sub["stain"].astype(str)) - {""}))
        testable = {t: v for t, v in groups.items() if len(v) >= 2}
        row = {**base, "n_treatments": len(groups), "replicates_per_treatment": ",".join(str(len(v)) for v in groups.values()),
               "welch_F": float("nan"), "welch_p": float("nan"), "note": ""}
        if len(testable) >= 2 and len(testable) == len(groups):
            try:
                row["welch_F"], row["welch_p"] = welch_anova(list(testable.values()))
            except Exception as exc:  # noqa: BLE001 - e.g. zero variance in every group
                row["note"] = f"not testable: {exc}"
        else:
            row["note"] = "needs >= 2 treatments with >= 2 replicates each"
        omni.append(row)
        ctrl = order[0] if order and is_control(order[0]) else ""
        pairs += _pair_rows({**base, "treatment": "", "compare": "treatment"}, groups, _select(long, GROUP, key),
                            "treatment", ctrl)

    # cell types (e.g. WT vs KO) within each staining panel and treatment
    cgroup = ["antibody", "metric", "metric_label", "feature", "treatment"]
    for key, sub in rep.groupby(cgroup, dropna=False):
        order = _ordered(set(sub["cell_type"]), orders["cell_type"])
        if len(order) < 2:
            continue
        groups = {c: sub.loc[sub["cell_type"] == c, "median"].dropna().to_numpy() for c in order}
        groups = {c: v for c, v in groups.items() if len(v)}
        base = {**dict(zip(cgroup, key)), "cell_type": "", "compare": "cell_type",
                "stain": ";".join(sorted(set(sub["stain"].astype(str)) - {""}))}
        ref = order[0] if is_wildtype(order[0]) else ""
        pairs += _pair_rows(base, groups, _select(long, cgroup, key), "cell_type", ref)
    omni_df = pd.DataFrame(omni, columns=OMNI_COLUMNS) if omni else pd.DataFrame(columns=OMNI_COLUMNS)
    omni_df["welch_q_bh"] = bh_qvalues(omni_df["welch_p"].to_numpy()) if len(omni_df) else []
    return omni_df, (pd.DataFrame(pairs, columns=PAIR_COLUMNS) if pairs else pd.DataFrame(columns=PAIR_COLUMNS))


def _p_text(omni: pd.DataFrame, ab, metric, feat, cell_types) -> str:
    """Welch ANOVA p across treatments, per cell type."""
    if not len(omni):
        return ""
    om = omni[(omni["antibody"] == ab) & (omni["metric"] == metric) & (omni["feature"] == feat)]
    parts = []
    for ct in cell_types:
        r = om[om["cell_type"] == ct]
        if len(r) and np.isfinite(r["welch_p"].iloc[0]):
            parts.append(f"{ct + ' ' if ct and len(cell_types) > 1 else ''}{r['welch_p'].iloc[0]:.3g}")
    return ("Welch p: " + " · ".join(parts)) if parts else ""


def _ylim(values: np.ndarray) -> tuple[float, float]:
    v = values[np.isfinite(values)]
    if not len(v):
        return 0.0, 1.0
    lo, hi = float(v.min()), float(v.max())
    if lo >= 0 and lo < 0.25 * hi:
        lo = 0.0
    pad = 0.04 * ((hi - lo) or 1.0)
    return lo - (pad if lo != 0 else 0.5 * pad), hi + pad


def _positions(t_order: list, ct_order: list, gap: float = 0.6) -> dict:
    """x position of each (treatment, cell type): cell types side by side within a treatment group."""
    n = max(1, len(ct_order))
    return {(t, c): i * (n + gap) + j for i, t in enumerate(t_order) for j, c in enumerate(ct_order or [""])}


def _group_axis(ax, t_order: list, ct_order: list, pos: dict) -> None:
    """Two-level x axis: cell type under each column, treatment underlined beneath its columns."""
    from matplotlib.transforms import offset_copy

    xs = list(pos.values())
    ax.set_xlim(min(xs, default=0) - 0.6, max(xs, default=0) + 0.6)
    if len(ct_order) <= 1:
        ax.set_xticks([pos[(t, ct_order[0] if ct_order else "")] for t in t_order])
        ax.set_xticklabels(t_order, fontsize=9)
        return
    ax.set_xticks([pos[(t, c)] for t in t_order for c in ct_order])
    ax.set_xticklabels([c for _ in t_order for c in ct_order], fontsize=8)
    line_tr = offset_copy(ax.get_xaxis_transform(), fig=ax.figure, y=-19, units="points")
    for t in t_order:
        a, b = pos[(t, ct_order[0])], pos[(t, ct_order[-1])]
        ax.plot([a - 0.4, b + 0.4], [0, 0], transform=line_tr, color="black", lw=0.8, clip_on=False)
        ax.annotate(t, xy=((a + b) / 2, 0), xycoords=("data", "axes fraction"), xytext=(0, -22),
                    textcoords="offset points", ha="center", va="top", fontsize=9)


def _box(ax) -> None:
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(1.0)
    ax.tick_params(axis="y", labelsize=8)


def _swarm(y: np.ndarray, x_per_pt: float, y_per_pt: float, dot_pt: float, half_width: float,
           rng: np.random.Generator) -> np.ndarray:
    """Beeswarm x offsets: values are binned by one dot height and spread sideways within a bin
    (centre out); the swarm is narrowed to `half_width` if it would overflow its column."""
    if not len(y):
        return np.zeros(0)
    bins = np.floor(y / (dot_pt * y_per_pt)).astype(np.int64)
    off = np.zeros(len(y))
    used: dict[int, int] = {}
    for idx in rng.permutation(len(y)):
        k = used.get(bins[idx], 0)
        used[bins[idx]] = k + 1
        off[idx] = ((k + 1) // 2) * (1 if k % 2 else -1) * dot_pt * x_per_pt
    m = np.abs(off).max()
    return off * (half_width / m) if m > half_width else off


def _replicate_keys(rep_sub: pd.DataFrame) -> tuple[dict, list, bool]:
    """Replicate -> legend key. Replicate labels shared across conditions (a Replicate column) keep their
    label; otherwise (one well per replicate) replicates are numbered 1, 2, 3 within each condition."""
    r = rep_sub.assign(_r=rep_sub["replicate"].astype(str))
    conds = r[["_r", "treatment", "cell_type"]].drop_duplicates()
    shared = bool(len(conds)) and bool(conds["_r"].duplicated().any())
    keys = {}
    if shared:
        for i, lab in enumerate(sorted(conds["_r"].unique())):
            keys.update({(t, c, lab): lab for t, c in conds.loc[conds["_r"] == lab, ["treatment", "cell_type"]].values})
        levels = sorted(conds["_r"].unique())
    else:
        for (t, c), g in conds.groupby(["treatment", "cell_type"]):
            for i, lab in enumerate(sorted(g["_r"])):
                keys[(t, c, lab)] = str(i + 1)
        levels = sorted(set(keys.values()), key=int)
    return keys, levels, shared


def _size(n: float, n_max: float) -> float:
    """Marker area (pt^2) of a replicate dot, proportional to its number of nuclei."""
    return max(18.0, 380.0 * n / max(n_max, 1))


def _legend_handles(levels: list, n_max: float) -> list:
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator

    hs = [Line2D([], [], ls="", marker="o", markersize=7, color=REPLICATE_COLORS[i % len(REPLICATE_COLORS)],
                 mec="black", mew=0.6, label=str(lv)) for i, lv in enumerate(levels)]
    ns = [v for v in MaxNLocator(3, integer=True).tick_values(0, n_max) if 0 < v <= n_max * 1.05] or [n_max]
    hs.append(Line2D([], [], ls="", label="$n$ (nuclei)"))
    hs += [Line2D([], [], ls="", marker="o", markersize=math.sqrt(_size(v, n_max)), mfc="white", mec="black",
                  mew=0.8, label=f"{int(v)}") for v in ns]
    return hs


def add_legend(target, handles: list, **kw) -> None:
    target.legend(handles=handles, title="Replicate", fontsize=8, title_fontsize=9, frameon=False,
                  labelspacing=1.1, handletextpad=0.8, borderaxespad=0.2, **kw)


def superplot(ax, long_sub: pd.DataFrame, rep_sub: pd.DataFrame, title: str, p_text: str = "",
              t_order: "list | None" = None, ct_order: "list | None" = None, legend: bool = True) -> list:
    """SuperPlot (Lord et al. 2020): every nucleus as a beeswarm dot colored by replicate; each replicate's
    median as a large dot of the same color, its area proportional to the number of nuclei; mean +/- SD
    of the replicate medians in black. Treatments are groups on the x axis, cell types columns within a
    group. Returns the legend handles."""
    t_order = [t for t in (t_order or level_order(long_sub, "treatment", first=is_control))
               if (rep_sub["treatment"] == t).any()]
    ct_order = [c for c in (ct_order or level_order(long_sub, "cell_type", first=is_wildtype))
                if (rep_sub["cell_type"] == c).any()]
    pos = _positions(t_order, ct_order)
    keys, levels, _ = _replicate_keys(rep_sub)
    color = {lv: REPLICATE_COLORS[i % len(REPLICATE_COLORS)] for i, lv in enumerate(levels)}
    n_max = float(rep_sub["n_nuclei"].max()) if len(rep_sub) else 1.0

    ax.set_ylim(*_ylim(long_sub["value"].to_numpy(dtype=float)))
    _group_axis(ax, t_order, ct_order, pos)
    fig = ax.figure
    bb = ax.get_position()
    w_pt, h_pt = bb.width * fig.get_figwidth() * 72, bb.height * fig.get_figheight() * 72
    (x0, x1), (y0, y1) = ax.get_xlim(), ax.get_ylim()
    x_per_pt, y_per_pt = (x1 - x0) / w_pt, (y1 - y0) / h_pt
    dot_s = 5.0  # nucleus dot area (pt^2)
    rng = np.random.default_rng(0)
    for (t, ct), x in pos.items():
        sel = long_sub[(long_sub["treatment"] == t) & (long_sub["cell_type"] == ct)].dropna(subset=["value"])
        rs = rep_sub[(rep_sub["treatment"] == t) & (rep_sub["cell_type"] == ct)].sort_values("replicate")
        if len(sel):
            y = sel["value"].to_numpy(dtype=float)
            cols = [color.get(keys.get((t, ct, str(r)), ""), "0.5") for r in sel["replicate"].astype(str)]
            off = _swarm(y - y0, x_per_pt, y_per_pt, math.sqrt(dot_s) * 1.05, 0.42, rng)
            order = rng.permutation(len(y))  # interleave replicates so none hides another
            ax.scatter((x + off)[order], y[order], s=dot_s, c=np.asarray(cols, dtype=object)[order], lw=0, zorder=2)
        meds = rs["median"].dropna()
        if len(meds) > 1:
            m, sd = meds.mean(), meds.std(ddof=1)
            ax.errorbar(x, m, yerr=sd, color="black", capsize=3, lw=1, zorder=3)
            ax.hlines(m, x - 0.15, x + 0.15, color="black", lw=1.2, zorder=3)
        for _, r in rs.sort_values("n_nuclei", ascending=False).iterrows():  # small dots on top
            ax.scatter([x], [float(r["median"])], s=_size(r["n_nuclei"], n_max), zorder=4,
                       color=color.get(keys.get((t, ct, str(r["replicate"])), ""), "0.5"), edgecolor="black", lw=0.8)
    ax.set_title(title + (f"\n{p_text}" if p_text else ""), fontsize=9)
    _box(ax)
    handles = _legend_handles(levels, n_max)
    if legend:
        add_legend(ax, handles, loc="upper left", bbox_to_anchor=(1.01, 1.0))
    return handles


def count_bins(high: float) -> list[tuple[float, float, str]]:
    """Spot-count classes scaled to the data (`high` ~ 95th percentile): (lo, hi, label)."""
    for edges in ([1, 2, 3, 4, 5, 6], [1, 3, 6, 11, 16], [1, 6, 11, 21, 31], [1, 11, 26, 51, 101]):
        if high < edges[-1] * 1.5:
            break
    out = [(0, 0, "0")]
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append((lo, hi - 1, str(lo) if lo == hi - 1 else f"{lo}-{hi - 1}"))
    out.append((edges[-1], math.inf, f"≥{edges[-1]}"))
    return out


def spot_count_distribution(long: pd.DataFrame) -> pd.DataFrame:
    """Fraction of QC-passed nuclei in each spot-count class, per staining panel, channel (spots) or channel
    pair (double-positive spots), cell type and treatment (nuclei pooled over replicates - descriptive)."""
    cols = ["antibody", "feature", "cell_type", "treatment", "spots", "n_nuclei", "fraction"]
    d = long[long["metric"].isin(["spots_n", "n_double"])].dropna(subset=["value"]) if len(long) else long
    rows = []
    for (ab, feat), sub in d.groupby(["antibody", "feature"], dropna=False):
        bins = count_bins(float(sub["value"].quantile(0.95)))
        for (ct, t), g in sub.groupby(["cell_type", "treatment"], dropna=False):
            v = g["value"].to_numpy()
            for lo, hi, label in bins:
                rows.append({"antibody": ab, "feature": feat, "cell_type": ct, "treatment": t, "spots": label,
                             "n_nuclei": len(v), "fraction": float(((v >= lo) & (v <= hi)).mean())})
    return pd.DataFrame(rows, columns=cols)


def _count_bars(ax, dist: pd.DataFrame, t_order: list, ct_order: list, title: str) -> None:
    """Stacked bars: fraction of nuclei with 0, 1-2, ... spots, on the same two-level axis as the SuperPlots."""
    import matplotlib.pyplot as plt

    labels = list(dict.fromkeys(dist["spots"]))
    cmap = plt.get_cmap("viridis")
    colors = {lb: cmap(i / max(1, len(labels) - 1)) for i, lb in enumerate(labels)}
    t_order = [t for t in t_order if (dist["treatment"] == t).any()]
    ct_order = [c for c in ct_order if (dist["cell_type"] == c).any()]
    pos = _positions(t_order, ct_order)
    first = next(iter(pos), None)
    for (t, ct), x in pos.items():
        g = dist[(dist["treatment"] == t) & (dist["cell_type"] == ct)]
        if not len(g):
            continue
        bottom = 0.0
        for lb in labels:
            f = float(g.loc[g["spots"] == lb, "fraction"].sum())
            ax.bar(x, f, 0.8, bottom=bottom, color=colors[lb], edgecolor="white", lw=0.5,
                   label=lb if (t, ct) == first else None)
            bottom += f
        ax.text(x, 1.01, f"n={int(g['n_nuclei'].iloc[0])}", ha="center", va="bottom", fontsize=6, color="0.3")
    ax.set_ylim(0, 1.07)
    _group_axis(ax, t_order, ct_order, pos)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("fraction of nuclei", fontsize=8)
    ax.set_title(title, fontsize=9)
    _box(ax)
    ax.legend(title="spots / nucleus", fontsize=8, title_fontsize=9, frameon=False, bbox_to_anchor=(1.01, 1.0),
              loc="upper left")


def write_superplots(long: pd.DataFrame, rep: pd.DataFrame, omni: pd.DataFrame, out_dir,
                     orders: "dict | None" = None) -> list:
    """Per staining panel, with all cell types in the same graphs: a spots-per-nucleus figure, an
    overview figure, and one PNG per metric. Also writes spot_count_distribution.csv next to
    `out_dir`. Returns the figure paths (spots-per-nucleus first)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..plate import slug

    orders = orders or default_orders(long)
    written = []
    out_dir.mkdir(parents=True, exist_ok=True)
    dist = spot_count_distribution(long)
    dist.to_csv(out_dir.parent / "spot_count_distribution.csv", index=False)
    sel = long[(long["metric"].isin(SUPERPLOT_CHANNEL) & (long["feature"].str.contains("×") == False))  # noqa: E712
               | (long["metric"].isin(SUPERPLOT_COLOC) & long["feature"].str.contains("×"))]
    for ab, lsub in sel.groupby("antibody", dropna=False):
        rsub = rep[rep["antibody"] == ab]
        t_order = _ordered(set(lsub["treatment"]), orders["treatment"])
        ct_order = _ordered(set(lsub["cell_type"]), orders["cell_type"])
        width = max(3.5, 0.8 + 0.42 * len(t_order) * (max(1, len(ct_order)) + 0.6))
        name = slug(ab) or "panel"
        cts = " vs ".join(c for c in ct_order if c)
        header = f"{ab}{' · ' + cts if cts else ''}"

        # spots per nucleus: SuperPlot + count distribution, one row per channel, then one per channel
        # pair for the double-positive spots
        rows = [(f, "spots_n") for f in dict.fromkeys(lsub.loc[lsub["metric"] == "spots_n", "feature"])]
        rows += [(f, "n_double") for f in dict.fromkeys(lsub.loc[lsub["metric"] == "n_double", "feature"])]
        if rows:
            fig, axes = plt.subplots(len(rows), 2, figsize=(2 * width + 3.2, 4.0 * len(rows)), squeeze=False)
            for (a1, a2), (feat, m) in zip(axes, rows):
                ls = lsub[(lsub["metric"] == m) & (lsub["feature"] == feat)]
                rs = rsub[(rsub["metric"] == m) & (rsub["feature"] == feat)]
                st = ";".join(sorted(set(ls["stain"].astype(str)) - {"", "×"})) or feat
                what = "spots per nucleus" if m == "spots_n" else "double-positive spots per nucleus"
                if m != "spots_n":
                    st, feat_txt = st.replace("×", " + "), feat.replace("×", "+")
                else:
                    feat_txt = feat
                superplot(a1, ls, rs, f"{st} ({feat_txt}): {what}", _p_text(omni, ab, m, feat, ct_order),
                          t_order, ct_order)
                a1.set_ylabel(what, fontsize=9)
                _count_bars(a2, dist[(dist["antibody"] == ab) & (dist["feature"] == feat)], t_order, ct_order,
                            f"{st} ({feat_txt}): distribution (nuclei pooled)")
            fig.suptitle(f"{header} - spots per nucleus", fontsize=11)
            fig.tight_layout()
            path = out_dir / f"{name}_spots_per_nucleus.png"
            fig.savefig(path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            written.append(path)

        keys = list(lsub.groupby(["metric", "feature"], sort=False).groups)
        keys.sort(key=lambda k: ("×" in k[1], k[1], (SUPERPLOT_CHANNEL + SUPERPLOT_COLOC).index(k[0])))
        if not keys:
            continue
        ncol = 4 if width <= 4.5 else 3 if width <= 6.5 else 2
        nrow = math.ceil(len(keys) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(width * ncol + 1.2, 3.9 * nrow), squeeze=False)
        for ax in axes.ravel()[len(keys):]:
            ax.axis("off")
        handles = []
        for ax, (m, feat) in zip(axes.ravel(), keys):
            ls = lsub[(lsub["metric"] == m) & (lsub["feature"] == feat)]
            rs = rsub[(rsub["metric"] == m) & (rsub["feature"] == feat)]
            st = ";".join(sorted(set(ls["stain"].astype(str)) - {"", "×"}))
            ptxt = _p_text(omni, ab, m, feat, ct_order)
            label = ls["metric_label"].iloc[0] if len(ls) else m
            title = f"{st} ({feat})\n{label}" if st else f"{feat}\n{label}"
            handles = superplot(ax, ls, rs, title, ptxt, t_order, ct_order, legend=False) or handles
            f1, a1 = plt.subplots(figsize=(width + 1.2, 3.9))
            f1.subplots_adjust(left=0.12, right=1 - 1.2 / (width + 1.2), bottom=0.2, top=0.85)
            superplot(a1, ls, rs, title, ptxt, t_order, ct_order)
            a1.set_ylabel(label, fontsize=9)
            p1 = out_dir / name / f"{slug(feat)}_{m}.png"
            p1.parent.mkdir(parents=True, exist_ok=True)
            f1.savefig(p1, dpi=200, bbox_inches="tight")
            plt.close(f1)
        fig.suptitle(f"{header} - dots: nuclei, colored by replicate; large dots: replicate medians "
                     "(area ∝ nuclei); black: mean ± SD of replicate medians", fontsize=11)
        fig.tight_layout(rect=(0, 0, 1 - 1.1 / (width * ncol + 1.2), 1))
        add_legend(fig, handles, loc="upper right", bbox_to_anchor=(1.0, 0.95))
        path = out_dir / f"{name}_overview.png"
        fig.savefig(path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        written.append(path)
    return written


def write_spot_brightness_plot(spots: pd.DataFrame, candidates: "pd.DataFrame | None", floors: dict,
                               channels: list[str], px_um: float, path) -> None:
    """Per channel: distribution of spot brightness (mean top-hat intensity, log scale) per treatment and
    cell type, with the fixed floor (solid) or, for per-nucleus channels, the floor an automatic fixed
    threshold would use (dashed). Two separated humps = dim grains vs real bodies."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from .segment import SEPARABILITY_MIN, auto_spot_threshold

    if not channels:
        return
    px2 = px_um * px_um
    fig, axes = plt.subplots(1, len(channels), figsize=(6 * len(channels), 4.2), squeeze=False)
    for ax, c in zip(axes[0], channels):
        if candidates is not None and (candidates["channel"] == c).any():
            by = candidates[candidates["channel"] == c]  # spots before the fixed floor
        else:
            d = spots[spots["channel"] == c] if len(spots) else spots
            if not len(d):
                ax.axis("off")
                continue
            by = d.assign(b=d["integrated_tophat"] / (d["area_um2"] / px2))
        by = by[by["b"] > 0]
        pooled = by["b"].to_numpy()
        if not pooled.size:
            ax.axis("off")
            continue
        bins = np.logspace(np.log10(pooled.min()), np.log10(pooled.max()), 50)
        ax.hist(pooled, bins=bins, color="0.85", label="all spots")
        styles = ["-", "--", ":", "-."]
        for i, t in enumerate(level_order(by, "treatment", first=is_control)):
            for j, ct in enumerate(level_order(by, "cell_type", first=is_wildtype)):
                v = by.loc[(by["treatment"] == t) & (by["cell_type"] == ct), "b"].to_numpy()
                if len(v):
                    h, e = np.histogram(v, bins=bins)
                    ax.step(e[:-1], h, where="post", color=REPLICATE_COLORS[i % len(REPLICATE_COLORS)],
                            ls=styles[j % len(styles)], lw=1.2, label=f"{t}{' · ' + ct if ct else ''}")
        if c in floors:
            f, sep = floors[c]["floor"], floors[c]["separability"]
            ax.axvline(f, color="black", lw=1.5, label=f"floor {f:.0f} ({floors[c]['mode']})")
            note = f"spots below {f:.0f} are dropped in every well"
        else:
            f, sep = auto_spot_threshold(pooled)
            if f:
                ax.axvline(f, color="black", lw=1.2, ls="--", label=f"automatic floor would be {f:.0f} (not applied)")
            note = ("per-nucleus threshold; one population - keep it" if not np.isfinite(sep) or sep < SEPARABILITY_MIN
                    else "per-nucleus threshold; two populations - consider Fixed")
        ax.set_xscale("log")
        ax.set_xlabel("spot brightness (mean top-hat intensity)", fontsize=9)
        ax.set_ylabel("spots", fontsize=9)
        sep_txt = f" · separability {sep:.2f}" if np.isfinite(sep) else ""
        ax.set_title(f"{c}: {note}{sep_txt}", fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        _box(ax)
    fig.suptitle("Spot brightness - two separated humps (dim grains vs bodies) support a fixed floor; one hump "
                 f"(separability below {SEPARABILITY_MIN}) does not", fontsize=10)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=130)
    plt.close(fig)
