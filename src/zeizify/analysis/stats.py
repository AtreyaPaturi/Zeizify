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
    "spot_mean_area_um2": "mean spot area (µm²)",
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
    "frac_a_overlap": "A spots overlapping B (fraction)",
    "frac_b_overlap": "B spots overlapping A (fraction)",
    "enrich_a": "A-in-B overlap vs chance",
    "nn_ab_um": "A→B nearest spot (µm)",
}
SUPERPLOT_CHANNEL = ["spots_n", "spots_area_frac", "cv", "enrichment"]
SUPERPLOT_COLOC = ["pearson", "tm1", "tm2", "frac_a_overlap"]
WILDTYPE_KEYWORDS = ("wt", "wildtype", "wild", "parental", "control", "ctrl")
CELL_TYPE_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3", "#937860", "#DA8BC3", "#8C8C8C"]
REPLICATE_MARKERS = "osD^vPX*<>"


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


def _xticks(ax, labels: list) -> None:
    long_labels = max((len(t) for t in labels), default=0) > 6
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=30 if long_labels else 0, ha="right" if long_labels else "center", fontsize=8)
    ax.set_xlim(-0.6, len(labels) - 0.4)


def superplot(ax, long_sub: pd.DataFrame, rep_sub: pd.DataFrame, title: str, p_text: str = "",
              t_order: "list | None" = None, ct_order: "list | None" = None, legend: bool = True) -> list:
    """Grouped SuperPlot (Lord et al. 2020): treatments on x, cell types side by side, one color per
    cell type. Nuclei are small jittered dots, replicate medians large dots, bar = mean +/- SD of the
    replicate medians. When replicate labels recur across conditions (biological replicates) each
    replicate gets its own marker. Returns the legend handles."""
    from matplotlib.lines import Line2D

    t_order = [t for t in (t_order or level_order(long_sub, "treatment", first=is_control))
               if (rep_sub["treatment"] == t).any()]
    ct_order = [c for c in (ct_order or level_order(long_sub, "cell_type", first=is_wildtype))
                if (rep_sub["cell_type"] == c).any()]
    reps = list(dict.fromkeys(rep_sub.sort_values("replicate")["replicate"].astype(str)))
    conds = rep_sub.assign(_r=rep_sub["replicate"].astype(str))[["_r", "treatment", "cell_type"]].drop_duplicates()
    shared = bool(len(conds)) and bool(conds["_r"].duplicated().any())
    markers = {r: (REPLICATE_MARKERS[i % len(REPLICATE_MARKERS)] if shared else "o") for i, r in enumerate(reps)}
    n = max(1, len(ct_order))
    slot = 0.8 / n
    rng = np.random.default_rng(0)
    for j, ct in enumerate(ct_order):
        color = CELL_TYPE_COLORS[j % len(CELL_TYPE_COLORS)]
        for i, t in enumerate(t_order):
            x = i + (j - (n - 1) / 2) * slot
            v = long_sub.loc[(long_sub["treatment"] == t) & (long_sub["cell_type"] == ct), "value"].dropna().to_numpy()
            if len(v):
                ax.scatter(x + rng.uniform(-0.3, 0.3, len(v)) * slot, v, s=4, color=color, alpha=0.25, lw=0)
            rs = rep_sub[(rep_sub["treatment"] == t) & (rep_sub["cell_type"] == ct)].sort_values("replicate")
            step = min(0.08, slot * 0.5 / max(1, len(rs)))
            for k, (_, r) in enumerate(rs.iterrows()):
                ax.scatter([x + (k - (len(rs) - 1) / 2) * step], [float(r["median"])], s=50, color=color,
                           marker=markers[str(r["replicate"])], edgecolor="black", lw=0.8, zorder=3)
            meds = rs["median"].dropna()
            if len(meds):
                m = meds.mean()
                ax.hlines(m, x - slot * 0.4, x + slot * 0.4, color="black", lw=1.5, zorder=4)
                if len(meds) > 1:
                    ax.errorbar(x, m, yerr=meds.std(ddof=1), color="black", capsize=3, lw=1, zorder=4)
    _xticks(ax, t_order)
    ax.set_title(title + (f"\n{p_text}" if p_text else ""), fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    handles = [Line2D([], [], ls="", marker="o", color=CELL_TYPE_COLORS[j % len(CELL_TYPE_COLORS)], mec="black",
                      label=ct or "cells") for j, ct in enumerate(ct_order)]
    if shared:
        handles += [Line2D([], [], ls="", marker=markers[r], color="white", mec="black", label=f"replicate {r}")
                    for r in reps]
    if legend and (len(ct_order) > 1 or shared):
        ax.legend(handles=handles, fontsize=7, frameon=False, loc="best")
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
    """Fraction of QC-passed nuclei in each spot-count class, per staining panel, channel, cell type and
    treatment (nuclei pooled over replicates - descriptive)."""
    cols = ["antibody", "feature", "cell_type", "treatment", "spots", "n_nuclei", "fraction"]
    d = long[long["metric"] == "spots_n"].dropna(subset=["value"]) if len(long) else long
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
    """Stacked bars: fraction of nuclei with 0, 1-2, ... spots; treatments on x, cell types side by side."""
    import matplotlib.pyplot as plt

    labels = list(dict.fromkeys(dist["spots"]))
    cmap = plt.get_cmap("viridis")
    colors = {lb: cmap(i / max(1, len(labels) - 1)) for i, lb in enumerate(labels)}
    t_order = [t for t in t_order if (dist["treatment"] == t).any()]
    ct_order = [c for c in ct_order if (dist["cell_type"] == c).any()]
    n = max(1, len(ct_order))
    slot = 0.8 / n
    for j, ct in enumerate(ct_order):
        for i, t in enumerate(t_order):
            g = dist[(dist["treatment"] == t) & (dist["cell_type"] == ct)]
            if not len(g):
                continue
            x = i + (j - (n - 1) / 2) * slot
            bottom = 0.0
            for lb in labels:
                f = float(g.loc[g["spots"] == lb, "fraction"].sum())
                ax.bar(x, f, slot * 0.9, bottom=bottom, color=colors[lb], edgecolor="white", lw=0.5,
                       label=lb if (i == 0 and j == 0) else None)
                bottom += f
            if n > 1:
                ax.text(x, 1.01, ct, ha="center", va="bottom", fontsize=7, color=CELL_TYPE_COLORS[j % len(CELL_TYPE_COLORS)],
                        fontweight="bold")
            ax.text(x, 0.5, f"n={int(g['n_nuclei'].iloc[0])}", ha="center", va="center", fontsize=6, color="white",
                    rotation=90)
    _xticks(ax, t_order)
    ax.set_ylim(0, 1.08 if n > 1 else 1.0)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("fraction of nuclei", fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.set_title(title, fontsize=9)
    ax.legend(title="spots / nucleus", fontsize=7, title_fontsize=7, frameon=False, bbox_to_anchor=(1.0, 1.0),
              loc="upper left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


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
        width = max(4.0, 1.1 + 0.45 * len(t_order) * len(ct_order))
        name = slug(ab) or "panel"
        cts = " vs ".join(c for c in ct_order if c)
        header = f"{ab}{' · ' + cts if cts else ''}"

        # spots per nucleus: SuperPlot + count distribution, one row per channel
        feats = list(dict.fromkeys(lsub.loc[lsub["metric"] == "spots_n", "feature"]))
        if feats:
            fig, axes = plt.subplots(len(feats), 2, figsize=(2 * width + 1.5, 3.9 * len(feats)), squeeze=False,
                                     gridspec_kw={"width_ratios": [1, 1.15]})
            for (a1, a2), feat in zip(axes, feats):
                ls = lsub[(lsub["metric"] == "spots_n") & (lsub["feature"] == feat)]
                rs = rsub[(rsub["metric"] == "spots_n") & (rsub["feature"] == feat)]
                st = ";".join(sorted(set(ls["stain"].astype(str)) - {""})) or feat
                superplot(a1, ls, rs, f"{st} ({feat}): spots per nucleus",
                          _p_text(omni, ab, "spots_n", feat, ct_order), t_order, ct_order)
                a1.set_ylabel("spots per nucleus", fontsize=8)
                a1.set_ylim(bottom=0)
                _count_bars(a2, dist[(dist["antibody"] == ab) & (dist["feature"] == feat)], t_order, ct_order,
                            f"{st} ({feat}): distribution (nuclei pooled)")
            fig.suptitle(f"{header} - spots per nucleus (dots: nuclei; large dots: replicate medians; "
                         "bar: mean ± SD)", fontsize=11)
            fig.tight_layout()
            path = out_dir / f"{name}_spots_per_nucleus.png"
            fig.savefig(path, dpi=150)
            plt.close(fig)
            written.append(path)

        keys = list(lsub.groupby(["metric", "feature"], sort=False).groups)
        keys.sort(key=lambda k: ("×" in k[1], k[1], (SUPERPLOT_CHANNEL + SUPERPLOT_COLOC).index(k[0])))
        if not keys:
            continue
        ncol = 4 if width <= 4.5 else 3 if width <= 6 else 2
        nrow = math.ceil(len(keys) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(width * ncol, 3.6 * nrow), squeeze=False)
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
            f1, a1 = plt.subplots(figsize=(width, 3.8))
            superplot(a1, ls, rs, title, ptxt, t_order, ct_order)
            f1.tight_layout()
            p1 = out_dir / name / f"{slug(feat)}_{m}.png"
            p1.parent.mkdir(parents=True, exist_ok=True)
            f1.savefig(p1, dpi=200)
            plt.close(f1)
        if len(handles) > 1:
            fig.legend(handles=handles, loc="upper right", fontsize=9, frameon=False)
        fig.suptitle(f"{header} - dots: nuclei; large dots: replicate medians; bar: mean ± SD", fontsize=11)
        fig.tight_layout(rect=(0, 0, 0.94 if len(handles) > 1 else 1, 1))
        path = out_dir / f"{name}_overview.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)
    return written
