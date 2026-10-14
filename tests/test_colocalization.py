"""Multi-stain (co-localization) plate maps: per-channel antibodies, panels and colors."""

import csv

import numpy as np
from PIL import Image

from conftest import FakeCzi, make_entries
from zeizify.contrast import default_control_choice, write_intensity_report
from zeizify.extraction import match_stain, process_all
from zeizify.layout import (
    candidate_panels,
    channel_color,
    compose_rgb,
    generate_png_layouts,
    panel_caption,
    parse_panel,
    stains_by_channel,
)
from zeizify.plate import load_plate_map, write_plate_map

CH = ["AF488", "AF594", "DAPI"]


def _coloc_map(tmp_path):
    p = tmp_path / "map.csv"
    write_plate_map(p, {
        "A1": {"CellType": "PAR", "Treatment": "Mock", "Antibody_AF488": "CD7", "Antibody_AF594": "HLA"},
        "A2": {"CellType": "PAR", "Treatment": "IFN", "Antibody_AF488": "CD7", "Antibody_AF594": "HLA"},
        "B1": {"CellType": "PAR", "Treatment": "Mock", "Antibody_AF488": "CD7"},  # single stain
    }, fluorophores=["AF488", "AF594"])
    return p


def test_csv_has_one_column_per_fluorophore(tmp_path):
    p = _coloc_map(tmp_path)
    header = p.read_text(encoding="utf-8").splitlines()[0]
    assert header == "Well,Row,Column,CellType,Treatment,Antibody_AF488,Antibody_AF594,Antibody,Label"
    m = load_plate_map(p)
    assert m["A1"]["Antibody"] == "CD7+HLA"
    assert m["A1"]["Label"] == "PAR_Mock_CD7_HLA"
    assert m["A1"]["Stains"] == {"AF488": "CD7", "AF594": "HLA"}
    assert m["B1"]["Antibody"] == "CD7" and m["B1"]["Stains"] == {"AF488": "CD7"}


def test_hand_made_csv_without_label_or_panel(tmp_path):
    p = tmp_path / "map.csv"
    p.write_text("Well,Treatment,Antibody_AF488,Antibody_AF647\nC3,DMSO,CD7,LAMP1\nC4,,,\n", encoding="utf-8")
    m = load_plate_map(p)
    assert list(m) == ["C3"]
    assert m["C3"]["Antibody"] == "CD7+LAMP1" and m["C3"]["Label"] == "DMSO_CD7_LAMP1"


def test_match_stain():
    stains = {"AF488": "CD7", "Alexa 594": "HLA"}
    assert match_stain("AF488", stains) == "CD7"
    assert match_stain("af-488", stains) == "CD7"
    assert match_stain("AF594", stains) == "HLA"  # by wavelength
    assert match_stain("DAPI", stains) == ""
    assert match_stain("AF488", {}) == ""


def test_extraction_tags_channels_and_groups_by_panel(tmp_path):
    out = tmp_path / "out"
    czi = FakeCzi(make_entries({0: (0, 0), 1: (0, 9000), 2: (9000, 0)}))
    pm = load_plate_map(_coloc_map(tmp_path))
    stats, skipped = process_all(czi, pm, {"A1": [0], "A2": [1], "B1": [2]}, CH, {0: "P1", 1: "P1", 2: "P1"}, out,
                                 log=lambda s: None)
    assert skipped == []
    tag = {(s.well, s.channel_name): s.stain for s in stats}
    assert tag[("A1", "AF488")] == "CD7" and tag[("A1", "AF594")] == "HLA" and tag[("A1", "DAPI")] == ""
    assert tag[("B1", "AF594")] == ""
    assert (out / "CD7_HLA" / "PAR_Mock_CD7_HLA_A1_P1_AF594.tif").exists()
    assert (out / "CD7" / "PAR_Mock_CD7_B1_P1_AF488.tif").exists()
    assert stains_by_channel(stats)["CD7+HLA"] == {"AF488": "CD7", "AF594": "HLA"}

    control = default_control_choice(stats)
    assert control["CD7+HLA"] == ("A1", "P1")
    rows = list(csv.DictReader(open(write_intensity_report(stats, control, out), encoding="utf-8")))
    assert {r["Stain"] for r in rows if r["Well"] == "A2" and r["Channel"] == "AF594"} == {"HLA"}

    panels = [("AF488", "AF594")]
    assert generate_png_layouts(stats, control, out, panels=panels, colors={"AF488": "Green", "AF594": "Magenta"}) == 3
    img = Image.open(out / "CD7_HLA" / "png" / "PAR_Mock_CD7_HLA_A1_P1_layout.png")
    assert img.size == (10, 8)  # exactly the one requested merge panel


def test_candidate_panels_cover_all_combinations():
    c = candidate_panels(["AF488", "AF594", "AF647", "DAPI"])
    assert ("AF488", "AF594") in c and ("AF594", "AF647") in c and ("AF488", "AF594", "AF647", "DAPI") in c
    assert len(c) == len(set(c))
    assert parse_panel("AF488 + AF594") == ("AF488", "AF594")


def test_colocalized_green_magenta_is_white():
    scaled = {"AF488": np.full((2, 2), 255, np.uint8), "AF594": np.full((2, 2), 255, np.uint8)}
    rgb = compose_rgb(scaled, ("AF488", "AF594"), {"AF488": "Green", "AF594": "Magenta"})
    assert rgb[0, 0].tolist() == [255, 255, 255]
    assert channel_color("AF594", {"AF594": "Cyan"}) == (0.0, 1.0, 1.0)
    assert panel_caption(("AF488",), {"AF488": "CD7"}) == "CD7 (AF488)"
    assert panel_caption(("AF488", "AF594", "DAPI"), {"AF488": "CD7", "AF594": "HLA"}) == "CD7 + HLA + DAPI"
