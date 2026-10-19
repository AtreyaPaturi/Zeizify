"""Zeizify Analyze: segmentation, spot detection and colocalization against synthetic ground truth;
statistics against reference implementations (pingouin, SciPy, scikit-image)."""

import numpy as np
import pandas as pd
import pytest
from skimage.measure import manders_coloc_coeff, pearson_corr_coeff

from analysis_synth import PX, make_field, write_zeizify_output
from zeizify.analysis.coloc import costes_thresholds, object_coloc, pixel_coloc, regression_tls
from zeizify.analysis.params import AnalysisParams
from zeizify.analysis.segment import detect_spots, segment_nuclei
from zeizify.analysis.stats import bh_qvalues, cliffs_delta, games_howell, hedges_g, welch_anova

P = AnalysisParams(nucleus_method="classical", spot_min_diameter_um=0.5, spot_max_diameter_um=3.0)


# ---------------- segmentation ----------------

def _grid_order(lab):
    from skimage.measure import regionprops

    props = sorted(regionprops(lab), key=lambda r: (round(r.centroid[0] / 100), r.centroid[1]))
    return [r.label for r in props]


def test_classical_nuclei_and_spot_counts_match_ground_truth():
    dapi, a, b, truth = make_field(seed=1)
    nuc = segment_nuclei(dapi, P, PX, "classical")
    assert nuc.max() == 9
    sa, _ = detect_spots(a, nuc, P, PX)
    sb, _ = detect_spots(b, nuc, P, PX)
    for lbl, (na, _, nb) in zip(_grid_order(nuc), truth):
        got_a = len(np.unique(sa[(nuc == lbl) & (sa > 0)]))
        got_b = len(np.unique(sb[(nuc == lbl) & (sb > 0)]))
        assert got_a == na and got_b == nb


def test_no_spots_called_in_a_diffuse_nucleus():
    dapi, a, _, _ = make_field(seed=2, n_spots=(0,) * 9, b_extra=0)
    nuc = segment_nuclei(dapi, P, PX, "classical")
    sa, _ = detect_spots(a, nuc, P, PX)
    assert sa.max() == 0  # the MAD noise floor keeps noise from being split into spots


def test_size_window_and_border_exclusion():
    dapi, *_ = make_field(seed=3)
    big = AnalysisParams(nucleus_method="classical", nucleus_diameter_um=40)  # expects far larger nuclei
    assert segment_nuclei(dapi, big, PX, "classical").max() == 0
    cut = np.roll(dapi, -100, axis=1)  # the left column of nuclei now straddles the image edge
    assert segment_nuclei(cut, P, PX, "classical").max() == 6
    keep = AnalysisParams(nucleus_method="classical", exclude_border_nuclei=False)
    assert segment_nuclei(cut, keep, PX, "classical").max() == 9


# ---------------- colocalization: reference implementations ----------------

def test_pearson_and_manders_match_scikit_image():
    rng = np.random.default_rng(0)
    a = rng.gamma(2, 50, 5000)
    b = 0.6 * a + rng.gamma(2, 30, 5000)
    r = pixel_coloc(a, b)
    assert r["pearson"] == pytest.approx(pearson_corr_coeff(a, b)[0], abs=1e-9)
    assert r["m1"] == pytest.approx(manders_coloc_coeff(a, b > 0), abs=1e-9)
    assert r["tm1"] == pytest.approx(manders_coloc_coeff(a, b >= r["costes_tb"]), abs=1e-9)
    assert r["tm2"] == pytest.approx(manders_coloc_coeff(b, a >= r["costes_ta"]), abs=1e-9)


def test_costes_threshold_zeroes_below_threshold_correlation():
    rng = np.random.default_rng(1)
    a = rng.gamma(2, 50, 20000)
    b = 0.8 * a + rng.normal(0, 20, 20000)
    ta, tb = costes_thresholds(a, b)
    m, c = regression_tls(a, b)
    assert tb == pytest.approx(m * ta + c)
    below = (a < ta) | (b < tb)
    assert abs(np.corrcoef(a[below], b[below])[0, 1]) < 0.02
    assert np.isnan(costes_thresholds(a, -a)[0])  # anti-correlated: undefined


def test_tls_regression_recovers_slope():
    rng = np.random.default_rng(2)
    a = rng.normal(100, 20, 50000)
    b = 2 * a + 5 + rng.normal(0, 1, 50000)
    m, c = regression_tls(a, b)
    assert m == pytest.approx(2, rel=0.01)


def test_icq_and_spearman_bounds():
    rng = np.random.default_rng(3)
    a, b = rng.random(3000), rng.random(3000)
    r = pixel_coloc(a, b)
    assert abs(r["icq"]) < 0.05 and abs(r["spearman"]) < 0.06
    r = pixel_coloc(a, a)
    assert r["icq"] == pytest.approx(0.5, abs=1e-6) and r["spearman"] == pytest.approx(1)


def test_object_coloc_counts_shared_spots():
    dapi, a, b, truth = make_field(seed=4, shared=0.5, b_extra=2)
    nuc = segment_nuclei(dapi, P, PX, "classical")
    sa, _ = detect_spots(a, nuc, P, PX)
    sb, _ = detect_spots(b, nuc, P, PX)
    for lbl, (na, ns, nb) in zip(_grid_order(nuc), truth):
        m = nuc == lbl
        oc = object_coloc(np.where(m, sa, 0), np.where(m, sb, 0), m, PX)
        assert oc["frac_a_overlap"] == pytest.approx(ns / na)
        assert oc["frac_b_overlap"] == pytest.approx(ns / nb)
        assert oc["n_a_overlap"] == oc["n_b_overlap"] == oc["n_double"] == ns
        assert (oc["n_a_only"], oc["n_b_only"]) == (na - ns, nb - ns)
        assert oc["enrich_a"] > 3  # shared spots overlap far more than chance
        assert oc["nn_ab_um"] < 0.5 or ns <= na / 2


def test_colocalized_vs_independent_pixels():
    dapi, a, b, _ = make_field(seed=5, shared=1.0, b_extra=0)
    _, a2, b2, _ = make_field(seed=6, n_spots=(5,) * 9, shared=0.0, b_extra=5)
    nuc = segment_nuclei(dapi, P, PX, "classical")
    m = nuc == 1
    assert pixel_coloc(a[m], b[m])["pearson"] > 0.9
    assert pixel_coloc(a2[m], b2[m])["pearson"] < 0.5


# ---------------- statistics: reference values from pingouin 0.7.0 ----------------

DATA = {"ctrl": np.array([10.1, 11.3, 9.8, 10.6]), "drugA": np.array([14.2, 16.8, 15.1]),
        "drugB": np.array([10.9, 13.5, 12.2, 11.1, 12.8])}


def test_games_howell_matches_pingouin():
    ref = {("ctrl", "drugA"): (5.924937, 2.745762, 0.024319), ("ctrl", "drugB"): (2.779101, 6.588931, 0.066070),
           ("drugA", "drugB"): (-3.594062, 3.712019, 0.054279)}
    for row in games_howell(DATA):
        t, df, p = ref[(row["group_a"], row["group_b"])]
        assert row["t"] == pytest.approx(t, abs=1e-5)
        assert row["df"] == pytest.approx(df, abs=1e-5)
        assert row["p"] == pytest.approx(p, abs=1e-4)


def test_welch_anova_matches_pingouin():
    F, p = welch_anova(list(DATA.values()))
    assert F == pytest.approx(16.360785, abs=1e-5)
    assert p == pytest.approx(0.007919, abs=1e-5)


def test_effect_sizes_and_bh():
    assert cliffs_delta(np.array([1.0, 2, 3]), np.array([4.0, 5, 6])) == 1
    assert cliffs_delta(np.array([1.0, 2, 3]), np.array([1.0, 2, 3])) == 0
    assert hedges_g(np.array([1.0, 2, 3]), np.array([1.0, 2, 3])) == 0
    from scipy.stats import false_discovery_control

    p = np.array([0.01, 0.04, 0.03, np.nan, 0.5, 0.002])
    q = bh_qvalues(p)
    ok = np.isfinite(p)
    assert np.allclose(q[ok], false_discovery_control(p[ok], method="bh")) and np.isnan(q[3])


# ---------------- end to end ----------------

def test_run_analysis_end_to_end(tmp_path):
    from zeizify.analysis.run import run_analysis

    design = []
    for i, (trt, shared) in enumerate([("Mock", 0.9), ("Mock", 0.9), ("Drug", 0.2), ("Drug", 0.2)]):
        for f in ("P1", "P2"):
            design.append({"well": f"A{i + 1}", "field": f, "treatment": trt, "seed": 10 * i + int(f[1]), "shared": shared})
    out = tmp_path / "out"
    write_zeizify_output(out, design)
    params = AnalysisParams(nucleus_method="classical", spot_min_diameter_um=0.5, workers=2)
    res = run_analysis(out, params, log=lambda s: None)
    adir = res.analysis_dir
    for name in ("nuclei", "spots", "colocalization", "fields", "replicate_summary", "condition_summary",
                 "stats_omnibus", "stats_pairwise"):
        assert (adir / f"{name}.csv").exists()
    assert (adir / "analysis_manifest.json").exists() and res.overview_plots and res.overview_plots[0].exists()
    assert res.overview_plots[1].name.endswith("_single_double_positive.png") and res.overview_plots[1].exists()
    assert res.overview_plots[2].name.endswith("_spot_categories.png") and res.overview_plots[2].exists()
    assert res.overview_plots[3].name.endswith("_spot_area_per_nucleus.png") and res.overview_plots[3].exists()
    assert (adir / "spot_categories.csv").exists()
    assert res.overview_plots[0].name.endswith("_spots_per_nucleus.png") and (adir / "spot_count_distribution.csv").exists()
    nuclei = pd.read_csv(adir / "nuclei.csv")
    assert len(nuclei) == 8 * 9 and nuclei["qc_pass"].all()
    assert set(nuclei["AF488_stain"]) == {"X"}
    omni = pd.read_csv(adir / "stats_omnibus.csv")
    row = omni[(omni["metric"] == "pearson")].iloc[0]
    assert row["replicates_per_treatment"] == "2,2" and row["welch_p"] < 0.05  # 90 % vs 20 % shared spots
    zero_var = omni[(omni["metric"] == "frac_a_overlap")].iloc[0]  # identical replicates: test undefined, not crash
    assert zero_var["replicates_per_treatment"] == "2,2"
    pairs = pd.read_csv(adir / "stats_pairwise.csv")
    pr = pairs[pairs["metric"] == "frac_a_overlap"].iloc[0]
    assert pr["group_a"] == "Mock" and pr["vs_control"] and pr["cliffs_delta_nuclei"] < -0.9
    # masks are cached: a second run with different spot settings does not re-segment nuclei
    lines = []
    run_analysis(out, AnalysisParams(nucleus_method="classical", spot_min_diameter_um=0.5, spot_sensitivity=8,
                                     workers=2, qc_overlays=False), log=lines.append)
    assert any("8 cached, 0 to segment" in ln for ln in lines)


def test_mitotic_like_nucleus_is_flagged(tmp_path):
    from zeizify.analysis.io import load_dataset
    from zeizify.analysis.measure import measure_field

    out = tmp_path / "o"
    write_zeizify_output(out, [{"well": "A1", "field": "P1", "treatment": "Mock", "seed": 1, "shared": 0.5}])
    ds = load_dataset(out)
    dapi, a, b, _ = make_field(seed=1, dapi_scale=[1, 1, 1, 1, 3, 1, 1, 1, 1])
    import tifffile

    tifffile.imwrite(ds.fields[0].channels["DAPI"], dapi)
    nuc = segment_nuclei(dapi, P, PX, "classical")
    r = measure_field(ds.fields[0], nuc, "DAPI", ["AF488", "AF594"], P, PX)
    flagged = [n for n in r.nuclei if not n["qc_pass"]]
    assert len(flagged) == 1 and "bright DAPI" in flagged[0]["qc_reason"]


def test_treatment_and_cell_type_order_follow_plate_map():
    from zeizify.analysis.stats import default_orders

    long = pd.DataFrame({"well": ["C2", "A2", "D2", "B2", "G2", "E2"],
                         "treatment": ["MG132", "MOCK", "ML792+MG132", "ML792", "MG132", "MOCK"],
                         "cell_type": ["WT", "WT", "WT", "WT", "KO1", "KO1"]})
    o = default_orders(long)
    assert o["treatment"] == ["MOCK", "ML792", "MG132", "ML792+MG132"]
    assert o["cell_type"] == ["WT", "KO1"]
    o = default_orders(long, treatments=["mg132"], cell_types=["KO1"])
    assert o["treatment"][0] == "MG132" and o["cell_type"] == ["KO1", "WT"]


def test_cell_types_compared_within_treatment_and_count_distribution():
    from zeizify.analysis.stats import replicate_summary, run_tests, spot_count_distribution

    rng = np.random.default_rng(1)
    rows = []
    for ct, shift in (("WT", 0), ("KO1", 5)):
        for t in ("Mock", "Drug"):
            for w in range(3):
                well = f"{'A' if ct == 'WT' else 'E'}{w + 1}{t}"
                for v in rng.poisson(10 + shift, 30):
                    rows.append({"antibody": "X", "cell_type": ct, "treatment": t, "replicate": well, "well": well,
                                 "metric": "spots_n", "metric_label": "spots per nucleus", "feature": "AF488",
                                 "stain": "X", "value": float(v)})
    long = pd.DataFrame(rows)
    _, pairs = run_tests(replicate_summary(long), long)
    ct_pairs = pairs[pairs["compare"] == "cell_type"]
    assert set(ct_pairs["treatment"]) == {"Mock", "Drug"}
    r = ct_pairs.iloc[0]
    assert (r["group_a"], r["group_b"]) == ("WT", "KO1") and r["vs_control"] and r["gh_p"] < 0.01
    dist = spot_count_distribution(long)
    sums = dist.groupby(["cell_type", "treatment"])["fraction"].sum()
    assert np.allclose(sums, 1.0) and set(dist["n_nuclei"]) == {90}


def test_double_positive_spot_counts_merge_touching_bodies():
    la, lb = np.zeros((40, 40), int), np.zeros((40, 40), int)
    la[5:15, 5:15] = 1  # one A body overlapped by two B foci -> 1 double-positive spot
    lb[5:9, 5:9], lb[11:15, 11:15] = 1, 2
    la[25:30, 25:30] = 2  # A and B overlapping once
    lb[27:32, 27:32] = 3
    la[2:4, 30:32] = 3  # A alone
    lb[34:37, 2:5] = 4  # B alone
    oc = object_coloc(la, lb, np.ones((40, 40), bool), 0.1)
    assert (oc["n_a"], oc["n_b"], oc["n_a_overlap"], oc["n_b_overlap"], oc["n_double"]) == (3, 4, 2, 3, 2)
    assert (oc["n_a_only"], oc["n_b_only"]) == (1, 1)


def test_fixed_spot_floor_drops_dim_spots():
    from zeizify.analysis.segment import spot_brightness

    dapi, a, b, truth = make_field(seed=5)
    nuc = segment_nuclei(dapi, P, PX, "classical")
    spots, th = detect_spots(a, nuc, P, PX)
    bright = spot_brightness(spots, th)
    floor = float(np.median(bright))
    kept, th2 = detect_spots(a, nuc, P, PX, min_brightness=floor)
    kb = spot_brightness(kept, th2)
    assert 0 < len(kb) < len(bright) and (kb >= floor).all()
    assert detect_spots(a, nuc, P, PX, min_brightness=1e12)[0].max() == 0


def test_halfmax_size_is_independent_of_brightness():
    from zeizify.analysis.segment import halfmax_areas

    yy, xx = np.mgrid[:60, :60]
    blob = np.exp(-((yy - 30) ** 2 + (xx - 30) ** 2) / (2 * 3.0 ** 2))
    nuc = np.ones((60, 60), int)
    sizes = {}
    for amp in (1.0, 4.0):
        th = amp * blob
        lab = (th > 0.2).astype(int)  # the same absolute threshold: a brighter spot gets a bigger outline
        sizes[amp] = (int(lab.sum()), halfmax_areas(lab, th, nuc, 10)[1])
    assert sizes[4.0][0] > 1.5 * sizes[1.0][0]
    assert sizes[4.0][1] == sizes[1.0][1] == pytest.approx(2 * np.log(2) * np.pi * 9, rel=0.1)


def test_pair_label_names_the_stains():
    from zeizify.analysis.stats import pair_label

    assert pair_label("A-only spots per nucleus", "PML×SUMO2/3") == "PML-only spots per nucleus"
    assert pair_label("A spots overlapping B (fraction)", "PML×SUMO2/3") == "PML spots overlapping SUMO2/3 (fraction)"
    assert pair_label("A→B nearest spot (µm)", "X×") == "A→B nearest spot (µm)"


def test_auto_spot_threshold_splits_two_populations_only():
    from zeizify.analysis.segment import SEPARABILITY_MIN, auto_spot_threshold

    rng = np.random.default_rng(0)
    grains, bodies = 10 ** rng.normal(2.8, 0.1, 3000), 10 ** rng.normal(3.5, 0.1, 2000)
    t, sep = auto_spot_threshold(np.concatenate([grains, bodies]))
    assert 10 ** 2.95 < t < 10 ** 3.35 and sep > 0.85
    _, sep1 = auto_spot_threshold(10 ** rng.normal(3.0, 0.2, 5000))
    assert sep1 < SEPARABILITY_MIN


def test_spot_category_table_shares_and_errors():
    from zeizify.analysis.stats import spot_category_table

    rows = []
    for rep_ in ("W1", "W2"):
        for m, vals in (("n_double", [4, 6]), ("n_a_only", [1, 3]), ("n_b_only", [0, 0])):
            for v in vals:
                rows.append({"antibody": "P", "cell_type": "WT", "treatment": "Mock", "replicate": rep_, "well": rep_,
                             "metric": m, "metric_label": m, "feature": "AF488×AF594", "stain": "PML×SUMO",
                             "value": float(v)})
    tab = spot_category_table(pd.DataFrame(rows)).set_index("category")
    assert tab.loc["n_double", "mean_per_nucleus"] == 5 and tab.loc["n_a_only", "mean_per_nucleus"] == 2
    assert tab.loc["n_double", "fraction_of_spots"] == pytest.approx(10 / 14)
    assert tab["fraction_of_spots"].sum() == pytest.approx(1.0)
    assert set(tab["error_type"]) == {"SD of replicate means"} and tab.loc["n_double", "error"] == 0
    assert tab.loc["n_a_only", "category_label"] == "PML only"
    one = spot_category_table(pd.DataFrame([r for r in rows if r["replicate"] == "W1"])).set_index("category")
    assert set(one["error_type"]) == {"SEM over nuclei"}
    assert one.loc["n_double", "error"] == pytest.approx(np.std([4, 6], ddof=1) / np.sqrt(2))
