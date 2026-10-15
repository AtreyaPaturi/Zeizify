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
Comparisons are made within each staining panel and cell type, so only identically stained wells
are compared.
"""

from __future__ import annotations

import math
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


def treatment_order(treatments) -> list[str]:
    """Control-like treatments first, then in order of appearance."""
    seen = list(dict.fromkeys(treatments))
    ctrl = [t for t in seen if any(k in t.lower() for k in CONTROL_KEYWORDS)]
    return ctrl + [t for t in seen if t not in ctrl]


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
PAIR_COLUMNS = GROUP + ["stain", "group_a", "group_b", "vs_control", "replicates_a", "replicates_b", "nuclei_a",
                        "nuclei_b", "median_nuclei_a", "median_nuclei_b", "diff_replicate_means",
                        "hedges_g_replicates", "cliffs_delta_nuclei", "gh_t", "gh_df", "gh_p", "note"]


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


def run_tests(rep: pd.DataFrame, long: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(omnibus table, pairwise table)."""
    omni, pairs = [], []
    for key, sub in rep.groupby(GROUP, dropna=False):
        order = treatment_order(sub["treatment"])
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

        # Pairwise: descriptive effect sizes always; Games-Howell p only with >= 2 replicates per group.
        sel = long
        for col, val in zip(GROUP, key):
            sel = sel[sel[col] == val]
        ctrl = order[0] if order and any(k in order[0].lower() for k in CONTROL_KEYWORDS) else ""
        gh = {(r["group_a"], r["group_b"]): r for r in games_howell(testable)} if len(testable) >= 2 else {}
        for a, b in combinations(list(groups), 2):
            va = sel.loc[sel["treatment"] == a, "value"].to_numpy(dtype=float)
            vb = sel.loc[sel["treatment"] == b, "value"].to_numpy(dtype=float)
            g = gh.get((a, b), {})
            pairs.append({**base, "group_a": a, "group_b": b, "vs_control": bool(ctrl) and ctrl in (a, b),
                          "replicates_a": len(groups[a]), "replicates_b": len(groups[b]),
                          "nuclei_a": int(np.isfinite(va).sum()), "nuclei_b": int(np.isfinite(vb).sum()),
                          "median_nuclei_a": float(np.nanmedian(va)) if len(va) else float("nan"),
                          "median_nuclei_b": float(np.nanmedian(vb)) if len(vb) else float("nan"),
                          "diff_replicate_means": float(groups[b].mean() - groups[a].mean()),
                          "hedges_g_replicates": hedges_g(groups[a], groups[b]),
                          "cliffs_delta_nuclei": cliffs_delta(va, vb),
                          "gh_t": g.get("t", float("nan")), "gh_df": g.get("df", float("nan")),
                          "gh_p": g.get("p", float("nan")),
                          "note": "" if g else "descriptive only: needs >= 2 replicates per treatment for a test"})
    omni_df = pd.DataFrame(omni, columns=OMNI_COLUMNS) if omni else pd.DataFrame(columns=OMNI_COLUMNS)
    omni_df["welch_q_bh"] = bh_qvalues(omni_df["welch_p"].to_numpy()) if len(omni_df) else []
    return omni_df, (pd.DataFrame(pairs) if pairs else pd.DataFrame(columns=PAIR_COLUMNS))


def superplot(ax, long_sub: pd.DataFrame, rep_sub: pd.DataFrame, title: str, p_text: str = "") -> None:
    """Nuclei as small jittered dots colored by replicate; replicate medians as large dots of the same
    color; bar = mean of replicate medians +/- SD (Lord et al. 2020)."""
    import matplotlib.pyplot as plt

    order = treatment_order(list(rep_sub["treatment"]) or list(long_sub["treatment"]))
    reps = list(dict.fromkeys(rep_sub.sort_values("replicate")["replicate"]))
    cmap = plt.get_cmap("tab10")
    colors = {r: cmap(i % 10) for i, r in enumerate(reps)}
    rng = np.random.default_rng(0)
    for i, t in enumerate(order):
        rs = [r for r in reps if ((rep_sub["treatment"] == t) & (rep_sub["replicate"] == r)).any()]
        for j, r in enumerate(rs):
            off = (j - (len(rs) - 1) / 2) * 0.12
            v = long_sub.loc[(long_sub["treatment"] == t) & (long_sub["replicate"] == r), "value"].dropna().to_numpy()
            if len(v):
                ax.scatter(i + off + rng.uniform(-0.05, 0.05, len(v)), v, s=4, color=colors[r], alpha=0.25, lw=0)
            med = rep_sub.loc[(rep_sub["treatment"] == t) & (rep_sub["replicate"] == r), "median"]
            if len(med):
                ax.scatter([i + off], [float(med.iloc[0])], s=60, color=colors[r], edgecolor="black", lw=0.8, zorder=3)
        meds = rep_sub.loc[rep_sub["treatment"] == t, "median"].dropna()
        if len(meds):
            m = meds.mean()
            ax.hlines(m, i - 0.3, i + 0.3, color="black", lw=1.5, zorder=4)
            if len(meds) > 1:
                ax.errorbar(i, m, yerr=meds.std(ddof=1), color="black", capsize=4, lw=1, zorder=4)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=30, ha="right", fontsize=8)
    ax.set_title(title + (f"\n{p_text}" if p_text else ""), fontsize=9)
    ax.tick_params(axis="y", labelsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def write_superplots(long: pd.DataFrame, rep: pd.DataFrame, omni: pd.DataFrame, out_dir) -> list:
    """One overview figure per staining panel and cell type, plus one PNG per metric."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from ..plate import slug

    written = []
    sel = long[(long["metric"].isin(SUPERPLOT_CHANNEL) & (long["feature"].str.contains("×") == False))  # noqa: E712
               | (long["metric"].isin(SUPERPLOT_COLOC) & long["feature"].str.contains("×"))]
    for (ab, ct), lsub in sel.groupby(["antibody", "cell_type"], dropna=False):
        keys = list(lsub.groupby(["metric", "feature"], sort=False).groups)
        keys.sort(key=lambda k: ("×" in k[1], k[1], (SUPERPLOT_CHANNEL + SUPERPLOT_COLOC).index(k[0])))
        if not keys:
            continue
        ncol = 4
        nrow = math.ceil(len(keys) / ncol)
        fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.4 * nrow), squeeze=False)
        for ax in axes.ravel()[len(keys):]:
            ax.axis("off")
        for ax, (m, feat) in zip(axes.ravel(), keys):
            ls = lsub[(lsub["metric"] == m) & (lsub["feature"] == feat)]
            rs = rep[(rep["antibody"] == ab) & (rep["cell_type"] == ct) & (rep["metric"] == m) & (rep["feature"] == feat)]
            st = ";".join(sorted(set(ls["stain"].astype(str)) - {"", "×"}))
            om = omni[(omni["antibody"] == ab) & (omni["cell_type"] == ct) & (omni["metric"] == m) & (omni["feature"] == feat)] if len(omni) else omni
            ptxt = f"Welch p = {om['welch_p'].iloc[0]:.3g}" if len(om) and np.isfinite(om["welch_p"].iloc[0]) else ""
            label = ls["metric_label"].iloc[0] if len(ls) else m
            title = f"{st or feat} ({feat})\n{label}" if st else f"{feat}\n{label}"
            superplot(ax, ls, rs, title, ptxt)
            f1, a1 = plt.subplots(figsize=(4, 3.6))
            superplot(a1, ls, rs, title, ptxt)
            f1.tight_layout()
            p1 = out_dir / slug(f"{ab}_{ct}") / f"{slug(feat)}_{m}.png"
            p1.parent.mkdir(parents=True, exist_ok=True)
            f1.savefig(p1, dpi=200)
            plt.close(f1)
        fig.suptitle(f"{ab}{' · ' + ct if ct else ''} - dots: nuclei; large dots: replicate medians; bar: mean ± SD",
                     fontsize=11)
        fig.tight_layout()
        path = out_dir / f"{slug(ab) or 'panel'}{'_' + slug(ct) if ct else ''}_overview.png"
        out_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150)
        plt.close(fig)
        written.append(path)
    return written
