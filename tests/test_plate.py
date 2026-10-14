import pytest

from zeizify.plate import PLATE_FORMATS, load_plate_map, make_label, normalize_well, slug, write_plate_map


def test_well_normalization():
    assert normalize_well("a01") == "A1"
    assert normalize_well(" H12 ") == "H12"
    with pytest.raises(ValueError):
        normalize_well("12A")


def test_plate_formats():
    assert len(PLATE_FORMATS[96].wells()) == 96
    assert PLATE_FORMATS[384].wells()[-1] == "P24"


def test_label_and_slug():
    assert make_label("KO 1", "IFN-γ 10ng", "anti/CD7") == "KO_1_IFN_10ng_anti_CD7"
    assert slug("../../evil") == "evil"


def test_round_trip(tmp_path):
    p = tmp_path / "map.csv"
    write_plate_map(p, {"A1": {"CellType": "PAR", "Treatment": "Mock", "Antibody": "CD7"},
                        "B3": {"CellType": "KO1", "Treatment": "", "Antibody": "CD7"}})
    m = load_plate_map(p)
    assert set(m) == {"A1", "B3"}
    assert m["A1"]["Label"] == "PAR_Mock_CD7"
    assert m["B3"]["Label"] == "KO1_CD7"


def test_hand_edited_label_cannot_escape_output_dir(tmp_path):
    p = tmp_path / "map.csv"
    p.write_text("Well,Label,Antibody\nA01,../../x,CD7\n", encoding="utf-8")
    assert load_plate_map(p) == {"A1": {"Well": "A1", "Label": "x", "Antibody": "CD7", "Treatment": "", "CellType": "",
                                        "Stains": {}}}


def test_duplicate_wells_rejected(tmp_path):
    p = tmp_path / "map.csv"
    p.write_text("Well,Label\nA1,a\nA01,b\n", encoding="utf-8")
    with pytest.raises(ValueError):
        load_plate_map(p)
