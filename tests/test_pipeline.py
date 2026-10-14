"""End-to-end CLI run against an in-memory fake CziFile (metadata + subblocks)."""

import csv
import json

import czifile
import pytest

from conftest import grid_positions, make_entries
from zeizify.cli import main
from zeizify.plate import write_plate_map


def _xml(truth: dict[int, str], with_array_name: bool) -> str:
    scenes = "".join(
        f'<Scene Index="{s}" Name="P{s % 3 + 1}">' + (f"<ArrayName>{w}</ArrayName>" if with_array_name else "") + "</Scene>"
        for s, w in truth.items()
    )
    return (
        "<ImageDocument><Metadata><Information><Image><ComponentBitCount>14</ComponentBitCount><Dimensions>"
        '<Channels><Channel Id="Channel:0" Name="AF488"/><Channel Id="Channel:1" Name="AF594"/>'
        '<Channel Id="Channel:2" Name="DAPI"/></Channels>'
        f"<S><Scenes>{scenes}</Scenes></S></Dimensions></Image></Information>"
        '<Scaling><Items><Distance Id="X"><Value>6.5e-07</Value></Distance></Items></Scaling>'
        "</Metadata></ImageDocument>"
    )


@pytest.fixture(params=[True, False], ids=["metadata", "geometry"])
def fake_czi(request, monkeypatch):
    positions, truth = grid_positions()
    entries = make_entries(positions, n_z=2, shape=(6, 6))
    xml = _xml(truth, request.param)

    class FakeCziFile:
        def __init__(self, path):
            self.filtered_subblock_directory = tuple(entries)

        def metadata(self):
            return xml

        def __enter__(self):
            return self

        def __exit__(self, *a):
            return False

    monkeypatch.setattr(czifile, "CziFile", FakeCziFile)
    return request.param, truth


def test_cli_extract_end_to_end(fake_czi, tmp_path):
    source_is_metadata, truth = fake_czi
    czi = tmp_path / "plate.czi"
    czi.write_bytes(b"")
    pm = tmp_path / "map.csv"
    write_plate_map(pm, {"A1": {"CellType": "PAR", "Treatment": "Mock", "Antibody": "CD7"},
                         "H12": {"CellType": "KO1", "Treatment": "IFN", "Antibody": "CD7"},
                         "D5": {"CellType": "PAR", "Treatment": "DMSO", "Antibody": "HLA"}})
    out = tmp_path / "out"
    assert main(["extract", str(czi), "-p", str(pm), "-o", str(out), "--yes",
                 "--range", "CD7:AF594=0:500"]) == 0

    rows = list(csv.DictReader(open(out / "_scene_mapping" / "scene_to_well.csv", encoding="utf-8")))
    assert {int(r["Scene"]): r["Well"] for r in rows} == truth
    assert rows[0]["Source"] == ("metadata" if source_is_metadata else "geometry")

    tifs = sorted(p.name for p in (out / "CD7").glob("*.tif"))
    assert len(tifs) == 2 * 3 * 3  # 2 wells x 3 fields x 3 channels
    assert "PAR_Mock_CD7_A1_P1_AF488.tif" in tifs
    assert len(list((out / "HLA" / "png").glob("*.png"))) == 3

    report = list(csv.DictReader(open(out / "intensity_report.csv", encoding="utf-8")))
    assert {r["ControlWell"] for r in report if r["Antibody"] == "HLA"} == {"D5"}
    assert {r["ControlMax"] for r in report if r["Antibody"] == "CD7" and r["Channel"] == "AF594"} == {"500"}

    manifest = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["mapping"]["n_wells"] == 96
    assert manifest["czi"]["pixel_size_um"] == pytest.approx(0.65)
    assert manifest["contrast"]["custom_ranges"] == [{"antibody": "CD7", "channel": "AF594", "min": 0, "max": 500}]
