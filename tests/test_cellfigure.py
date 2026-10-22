"""Single-cell figure: cell ranking, cropping, layout and export on a synthetic analysed plate."""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pytest
from PIL import Image

from analysis_synth import PX, write_zeizify_output
from zeizify import cellfigure as cf
from zeizify.analysis.params import AnalysisParams


@pytest.fixture(scope="module")
def analysed(tmp_path_factory):
    from zeizify.analysis.run import run_analysis

    out = tmp_path_factory.mktemp("fig") / "out"
    design, i = [], 0
    for ct, row in (("WT", "A"), ("KO", "B")):
        for t, shared, col in (("Mock", 0.9, 1), ("Drug", 0.2, 2)):
            for f in ("P1", "P2"):
                design.append({"well": f"{row}{col}", "field": f, "treatment": t, "cell_type": ct, "seed": i,
                               "shared": shared})
                i += 1
    write_zeizify_output(out, design)
    run_analysis(out, AnalysisParams(nucleus_method="classical", spot_min_diameter_um=0.5, workers=2,
                                     qc_overlays=False), log=lambda s: None)
    return out


def test_conditions_and_ranking(analysed):
    cd = cf.load_cells(analysed)
    assert cd.conditions() == [("WT", "Mock"), ("WT", "Drug"), ("KO", "Mock"), ("KO", "Drug")]
    assert cd.signal == ["AF488", "AF594"] and cd.nuclear == "DAPI"
    r = cf.rank_candidates(cd, "KO", "Drug", 20, 20, n=5)
    assert 0 < len(r) <= 5
    assert set(r["cell_type"]) == {"KO"} and set(r["treatment"]) == {"Drug"}
    # counts match the condition medians (all synthetic nuclei have 6 spots in A)
    assert (r["AF488_spots_n"] - r["median_AF488_spots_n"]).abs().max() <= 1
    assert r["relax"].is_monotonic_increasing


def test_render_layout_and_shared_ranges(analysed):
    cd = cf.load_cells(analysed)
    picks = {k: cf.choice_from_row(v.iloc[0]) for k, v in cf.auto_choices(cd, cd.conditions(), 20, 16, n=2).items()}
    opts = cf.FigureOptions(panels=cf.parse_order("DAPI | AF488 | AF594 | AF488+AF594"), cell_types=["WT", "KO"],
                            treatments=["Mock", "Drug"], crop_w_um=20, crop_h_um=16, width_in=7.0)
    fig, crops, rng = cf.render(cd, picks, opts)
    w, h = fig.get_size_inches()
    assert w == pytest.approx(7.0)
    crop = next(iter(crops.values()))
    assert crop["own"].shape == (round(16 / PX), round(20 / PX))  # rows x cols = height x width
    assert crop["own"].any()
    assert set(rng) == {"DAPI", "AF488", "AF594"}  # one range per channel, shared by all conditions
    n_images = sum(1 for ax in fig.axes if ax.images)
    assert n_images == 4 * 2 * 2  # panels x treatments x cell types
    opts.rows, opts.blocks = "panel", "stacked"
    fig2, _, _ = cf.render(cd, picks, opts)
    assert fig2.get_size_inches()[1] > h  # stacked is taller
    opts.contrast = "panel"
    _, _, rng2 = cf.render(cd, picks, opts)
    assert all(isinstance(k, tuple) for k in rng2)


def test_export_writes_figure_table_legend(analysed):
    cd = cf.load_cells(analysed)
    ranked = cf.auto_choices(cd, cd.conditions(), 20, 20, n=3)
    picks = {k: cf.choice_from_row(v.iloc[0]) for k, v in ranked.items()}
    alt = ranked[("KO", "Mock")].iloc[1]
    picks[("KO", "Mock")] = cf.choice_from_row(alt, "manual")
    opts = cf.FigureOptions(panels=cf.default_figure_panels(cd.ds.channels), cell_types=["WT", "KO"],
                            treatments=["Mock", "Drug"], treatment_labels={"Drug": "Drug X"}, formats=["png", "pdf", "svg"],
                            dpi=150, scale_bar_um=5)
    res = cf.export(cd, picks, opts, "test figure")
    folder = res["folder"]
    assert {p.suffix for p in res["files"]} == {".png", ".pdf", ".svg"} and all(p.exists() for p in res["files"])
    png = Image.open(next(p for p in res["files"] if p.suffix == ".png"))
    assert abs(png.width - opts.width_in * 150) <= 1 and png.mode == "RGB"
    assert round(png.info["dpi"][0]) == 150
    assert "Drug X" in (folder / f"{folder.name}.svg").read_text(encoding="utf-8")  # editable, renamed label
    cells = pd.read_csv(folder / "cells.csv")
    row = cells[(cells["cell_type"] == "KO") & (cells["treatment"] == "Mock")].iloc[0]
    assert row["chosen"] == "manual" and int(row["nucleus_id"]) == int(alt["nucleus_id"])
    legend = (folder / "legend.txt").read_text(encoding="utf-8")
    assert "Scale bar, 5 µm" in legend and "identical in all conditions" in legend
    settings = json.loads((folder / "figure_settings.json").read_text(encoding="utf-8"))
    assert settings["options"]["panels"] == ["DAPI", "AF488", "AF594", "AF488+AF594"]


def test_parsers_and_defaults():
    assert cf.parse_order("DAPI | AF488+AF594, AF594") == [("DAPI",), ("AF488", "AF594"), ("AF594",)]
    assert cf.parse_labels("ML=ML792; MG = MG132 ;bad") == {"ML": "ML792", "MG": "MG132"}
    assert cf.default_figure_panels(["AF488", "AF594", "DAPI"]) == [("DAPI",), ("AF488",), ("AF594",),
                                                                    ("AF488", "AF594")]


def test_crop_pads_outside_the_image():
    a = np.arange(16, dtype=np.float32).reshape(4, 4)
    c = cf._crop(a, -1, 2, 3, 3)
    assert c.shape == (3, 3) and c[0].tolist() == [0, 0, 0] and c[1, :2].tolist() == [2, 3] and c[1, 2] == 0
