"""Plate-map editor state: dynamic fluorophores and migration of older saved maps."""

import json

from zeizify.gui.state import PlateState
from zeizify.plate import load_plate_map


def test_new_state_has_two_fluorophores(tmp_path):
    s = PlateState(tmp_path / "s.json")
    assert s.fluorophores == ["AF488", "AF594"]
    assert s.categories == ["cellType", "treatment", "stain:AF488", "stain:AF594"]


def test_add_rename_remove_and_export(tmp_path):
    s = PlateState(tmp_path / "s.json")
    assert s.add_fluorophore("AF647") and not s.add_fluorophore("af647")
    s.cats["stain:AF488"]["byWell"]["A1"] = s.add_label("stain:AF488", "CD7")["id"]
    s.cats["stain:AF647"]["byWell"]["A1"] = s.add_label("stain:AF647", "LAMP1")["id"]
    assert s.rename_fluorophore("AF647", "Cy5")
    s.remove_fluorophore("AF594")
    s.save()
    p = tmp_path / "m.csv"
    p.write_text(s.csv_text(), encoding="utf-8")
    m = load_plate_map(p)
    assert m["A1"]["Stains"] == {"AF488": "CD7", "Cy5": "LAMP1"} and m["A1"]["Antibody"] == "CD7+LAMP1"
    again = PlateState(tmp_path / "s.json")
    assert again.fluorophores == ["AF488", "Cy5"] and again.name_at("stain:Cy5", "A1") == "LAMP1"
    assert s.palette("stain:AF488")[0] != s.palette("stain:Cy5")[0]  # channels start on different hues


def test_v1_single_antibody_state_migrates(tmp_path):
    f = tmp_path / "s.json"
    lbl = {"id": "x", "name": "CD7", "color": "#000000"}
    f.write_text(json.dumps({"plateWells": 96, "cellType": {"labels": [], "byWell": {}},
                             "treatment": {"labels": [], "byWell": {}},
                             "antibody": {"labels": [lbl], "byWell": {"B2": "x"}}}), encoding="utf-8")
    s = PlateState(f)
    assert s.fluorophores == ["AF594"]
    assert s.name_at("stain:AF594", "B2") == "CD7"


def test_import_old_and_new_csv(tmp_path):
    old = tmp_path / "old.csv"
    old.write_text("Well,CellType,Treatment,Antibody,Label\nA1,PAR,Mock,CD7,PAR_Mock_CD7\n", encoding="utf-8")
    s = PlateState(tmp_path / "s.json")
    s.import_csv(old)
    assert s.name_at("stain:AF488", "A1") == "CD7"  # first declared fluorophore receives the single antibody
    new = tmp_path / "new.csv"
    new.write_text("Well,Antibody_AF555,Antibody_AF647\nA1,CD7,HLA\n", encoding="utf-8")
    s.import_csv(new)
    assert s.fluorophores == ["AF555", "AF647"] and s.name_at("stain:AF647", "A1") == "HLA"
