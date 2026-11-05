"""Whole cells / cytoplasm: marker detection, cell-to-nucleus assignment, classical segmentation, measurement
and the run integration, on synthetic fields."""

from __future__ import annotations

import csv
import json

import numpy as np
import pandas as pd
import tifffile

from analysis_synth import PX, SHAPE, make_field, write_zeizify_output
from zeizify.analysis.cells import assign_cells, cell_marker_for, drop_border_cells, segment_cells
from zeizify.analysis.params import AnalysisParams
from zeizify.analysis.segment import segment_nuclei

CENTERS = [(110 + 210 * i, 110 + 210 * j) for i in range(3) for j in range(3)]  # as in make_field
CELL_R = 95  # px: cells 47.5 µm across, 20 px apart
P = AnalysisParams(nucleus_method="classical", cell_method="classical", spot_min_diameter_um=0.5)


def _cell_bodies() -> np.ndarray:
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    lab = np.zeros(SHAPE, np.int32)
    for k, (cy, cx) in enumerate(CENTERS, start=1):
        lab[(yy - cy) ** 2 + (xx - cx) ** 2 <= CELL_R ** 2] = k
    return lab


def _marker(lab: np.ndarray, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    img = np.where(lab > 0, 900.0, 100.0) + rng.normal(0, 15, SHAPE)
    return np.clip(img, 0, 65535).astype(np.uint16)


def test_marker_from_stain_names_or_forced():
    stains = {"AF488": "F-Actin", "AF594": "Cyclin A"}
    chans = ["AF488", "AF594", "DAPI"]
    assert cell_marker_for(stains, chans, AnalysisParams()) == "AF488"
    assert cell_marker_for({"AF488": "Phalloidin"}, chans, AnalysisParams()) == "AF488"
    assert cell_marker_for({"AF488": "PML", "AF594": "SUMO2/3"}, chans, AnalysisParams()) is None
    assert cell_marker_for(stains, chans, AnalysisParams(cell_segmentation="off")) is None
    assert cell_marker_for({}, chans, AnalysisParams(cell_marker="AF594")) == "AF594"
    assert cell_marker_for({}, chans, AnalysisParams(cell_marker="AF647")) is None


def test_assign_cells_splits_merged_cells_and_drops_empty_ones():
    nuclei = np.zeros((20, 40), np.int32)
    nuclei[8:12, 4:8] = 1
    nuclei[8:12, 30:34] = 2
    nuclei[8:12, 18:21] = 3  # a nucleus with no cell around it
    cells = np.zeros_like(nuclei)
    cells[2:18, 1:14] = 5  # one "cell" that actually holds nuclei 1 and 2 (two touching cells found as one) ...
    cells[2:18, 26:39] = 5
    cells[0:5, 16:24] = 7  # ... and an anucleate fragment
    out = assign_cells(cells, nuclei)
    assert set(np.unique(out)) == {0, 1, 2, 3}
    assert (out[2:18, 1:14] == 1).all() and (out[2:18, 26:39] == 2).all()  # split by nearest nucleus
    assert (out[0:5, 16:24][nuclei[0:5, 16:24] == 0] == 0).all()  # fragment dropped
    assert (out == 3).sum() == (nuclei == 3).sum()  # only its own pixels
    # a cell touching the image edge keeps only its nucleus
    touching = np.where(out == 1, 1, 0)
    touching[0, 1:14] = 1  # extend cell 1 to the top edge
    assert (drop_border_cells(touching, nuclei) == 1).sum() == (nuclei == 1).sum()
    assert (drop_border_cells(np.where(out == 2, 2, 0), nuclei) == 2).sum() > (nuclei == 2).sum()  # inside: kept


def test_classical_cells_match_the_bodies():
    dapi, a, b, _ = make_field(seed=3)
    nuclei = segment_nuclei(dapi, P, PX, "classical")
    truth = _cell_bodies()
    cells = segment_cells(_marker(truth), dapi, nuclei, P, PX, "classical")
    assert set(np.unique(cells)) - {0} == set(np.unique(nuclei)) - {0}  # one cell per nucleus
    for nid in np.unique(nuclei)[1:]:
        t = truth[nuclei == nid][0]
        mine, ref = cells == nid, truth == t
        iou = (mine & ref).sum() / (mine | ref).sum()
        assert iou > 0.85, (nid, iou)


def test_measure_cytoplasm_and_nc_ratio(tmp_path):
    from zeizify.analysis.io import load_dataset
    from zeizify.analysis.measure import measure_field

    out = tmp_path / "o"
    write_zeizify_output(out, [{"well": "A1", "field": "P1", "treatment": "Mock", "seed": 1, "shared": 0.5}])
    ds = load_dataset(out)
    fi = ds.fields[0]
    dapi = tifffile.imread(fi.channels["DAPI"])
    truth = _cell_bodies()
    b = tifffile.imread(fi.channels["AF594"]).astype(np.float64)
    b[(truth > 0)] += 100  # cytoplasm 100 above background; nuclei were 300 above (plus spots)
    tifffile.imwrite(fi.channels["AF594"], np.clip(b, 0, 65535).astype(np.uint16))
    tifffile.imwrite(fi.channels["AF488"], _marker(truth))
    nuclei = segment_nuclei(dapi, P, PX, "classical")
    cells = segment_cells(_marker(truth), dapi, nuclei, P, PX, "classical")
    r = measure_field(fi, nuclei, "DAPI", ["AF594"], P, PX, cells=cells, cell_marker="AF488")
    df = pd.DataFrame(r.nuclei)
    assert len(df) == 9 and (df["cyto_area_um2"] > 0).all()
    area = np.pi * CELL_R ** 2 * PX ** 2
    assert df["cell_area_um2"].between(0.85 * area, 1.15 * area).all()
    assert df["AF594_cyto_mean_bgsub"].between(80, 120).all()  # ~100 above the cell-free background
    assert (df["AF594_nc_ratio"] > 2.5).all()
    assert "AF488_cyto_mean_bgsub" in df  # the marker channel is measured too
    assert r.field_row["cells_with_cytoplasm"] == 9


def test_run_segments_cells_only_where_a_marker_is_stained(tmp_path):
    from zeizify.analysis.run import run_analysis

    out = tmp_path / "o"
    write_zeizify_output(out, [{"well": "A1", "field": f, "treatment": "Mock", "seed": i, "shared": 0.5}
                               for i, f in enumerate(("P1", "P2"))])
    # field P1: AF488 becomes an F-actin stain with cell bodies; P2 keeps its nuclear stain "X"
    rep = out / "intensity_report.csv"
    rows = list(csv.DictReader(open(rep, encoding="utf-8")))
    for r in rows:
        if r["Field"] == "P1" and r["Channel"] == "AF488":
            r["Stain"] = "F-Actin"
            tifffile.imwrite(r["Path"], _marker(_cell_bodies()))
    with open(rep, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    res = run_analysis(out, AnalysisParams(nucleus_method="classical", cell_method="classical", workers=2,
                                           spot_min_diameter_um=0.5, qc_overlays=True), log=lambda s: None)
    n = pd.read_csv(res.analysis_dir / "nuclei.csv")
    p1, p2 = n[n["field"] == "P1"], n[n["field"] == "P2"]
    assert (p1["cyto_area_um2"] > 0).all() and (p1["cell_marker"] == "AF488").all()
    assert p2["cyto_area_um2"].isna().all()  # no marker in P2: nuclei only
    man = json.loads((res.analysis_dir / "analysis_manifest.json").read_text())
    assert man["cell_fields"] == 1 and man["cell_markers"] == ["AF488"] and man["cell_method"] == "classical"
    assert (res.analysis_dir / "masks").glob("*_cells.tif")
    assert len(list((res.analysis_dir / "masks").glob("*_cells.tif"))) == 1
    long_metrics = set(pd.read_csv(res.analysis_dir / "replicate_summary.csv")["metric"])
    assert {"nc_ratio", "cyto_mean_bgsub", "cell_area_um2"} <= long_metrics
    # cached: a second run does not segment the cells again
    lines = []
    run_analysis(out, AnalysisParams(nucleus_method="classical", cell_method="classical", workers=2,
                                     spot_min_diameter_um=0.5, qc_overlays=False), log=lines.append)
    assert any("1 cached, 0 to segment" in ln for ln in lines if ln.startswith("Cells:"))


def test_cells_absent_when_off(tmp_path):
    out = tmp_path / "o"
    write_zeizify_output(out, [{"well": "A1", "field": "P1", "treatment": "Mock", "seed": 2, "shared": 0.5}])
    from zeizify.analysis.run import run_analysis

    res = run_analysis(out, AnalysisParams(nucleus_method="classical", cell_segmentation="off", cell_marker="AF488",
                                           workers=1, spot_min_diameter_um=0.5, qc_overlays=False), log=lambda s: None)
    assert "cyto_area_um2" not in pd.read_csv(res.analysis_dir / "nuclei.csv")
