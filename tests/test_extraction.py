import csv

import numpy as np
import pytest
import tifffile
from PIL import Image

from conftest import FakeCzi, FakeEntry, make_entries
from zeizify.contrast import build_control_range, default_control_choice, write_intensity_report
from zeizify.extraction import ExtractionError, decode_scene_channel_mips, process_all
from zeizify.layout import channel_color, default_panels, generate_png_layouts

CH = ["AF488", "AF594", "DAPI"]
PLATE_MAP = {
    "A1": {"Well": "A1", "Label": "PAR_Mock_CD7", "Antibody": "CD7", "Treatment": "Mock", "CellType": "PAR"},
    "A2": {"Well": "A2", "Label": "PAR_IFN_CD7", "Antibody": "CD7", "Treatment": "IFN", "CellType": "PAR"},
}
SCENES = {"A1": [0, 1], "A2": [2]}
FIELDS = {0: "P1", 1: "P2", 2: "P1"}


@pytest.fixture
def czi():
    return FakeCzi(make_entries({0: (0, 0), 1: (0, 800), 2: (0, 9000)}))


def test_mip_is_max_over_z_per_channel(czi):
    mips = decode_scene_channel_mips(czi, 2)
    assert sorted(mips) == [0, 1, 2]
    for c, plane in mips.items():
        assert plane.shape == (8, 10)
        assert plane[0, 0] == 5000 + c  # the bright z=1 pixel survives
        assert plane[3, 3] == 2 * 100 + c * 10 + 3  # elsewhere the top z wins


def test_mosaic_scene_is_refused():
    e = make_entries({0: (0, 0)}, n_c=1, n_z=1)
    e.append(FakeEntry(0, ("C", "Z", "Y", "X", "S"), (0, 0, 0, 512, 0), e[0].data))
    with pytest.raises(ExtractionError, match="tile"):
        decode_scene_channel_mips(FakeCzi(e), 0)


def test_process_all_writes_scaled_tiffs_and_resumes(czi, tmp_out):
    logs = []
    stats, skipped = process_all(czi, PLATE_MAP, SCENES, CH, FIELDS, tmp_out, bit_ceiling=5001, px_size_um=0.325,
                                 log=logs.append)
    assert skipped == []
    assert len(stats) == 3 * 3
    p = tmp_out / "CD7" / "PAR_Mock_CD7_A1_P2_AF594.tif"
    assert p.exists()
    with tifffile.TiffFile(p) as t:
        assert t.is_imagej
        xres = t.pages[0].tags["XResolution"].value
        assert abs(xres[0] / xres[1] - 1 / 0.325) < 1e-3
    sat = {(s.well, s.field, s.channel_name): s.saturated for s in stats}
    assert sat[("A1", "P1", "AF594")] and not sat[("A1", "P1", "AF488")]

    mtime = p.stat().st_mtime_ns
    stats2, _ = process_all(czi, PLATE_MAP, SCENES, CH, FIELDS, tmp_out, log=logs.append)
    assert p.stat().st_mtime_ns == mtime  # resumed, not rewritten
    assert [(s.well, s.field, s.channel_name, s.max) for s in stats2] == [(s.well, s.field, s.channel_name, s.max) for s in stats]


def test_well_without_scenes_is_reported(czi, tmp_out):
    pm = dict(PLATE_MAP, B1={"Well": "B1", "Label": "x", "Antibody": "CD7"})
    _, skipped = process_all(czi, pm, SCENES, CH, FIELDS, tmp_out, log=lambda s: None)
    assert skipped == [("B1", "no scenes mapped to this well")]


def test_contrast_report_and_layouts(czi, tmp_out):
    stats, _ = process_all(czi, PLATE_MAP, SCENES, CH, FIELDS, tmp_out, log=lambda s: None)
    control = default_control_choice(stats)
    assert control == {"CD7": ("A1", "P1")}  # "Mock" is recognised as the control
    custom = {("CD7", "AF594"): (0, 100)}
    rng = build_control_range(stats, control, custom)
    assert rng[("CD7", "AF594")] == (0, 100)
    assert rng[("CD7", "DAPI")] == (23, 5002)  # MIP min is the top z value

    report = write_intensity_report(stats, control, tmp_out, custom)
    rows = list(csv.DictReader(open(report, encoding="utf-8")))
    assert len(rows) == 9
    af594 = [r for r in rows if r["Channel"] == "AF594"]
    assert all(r["ControlWell"] == "custom" and r["ControlMax"] == "100" for r in af594)

    n = generate_png_layouts(stats, control, tmp_out, custom)
    assert n == 3
    img = Image.open(tmp_out / "CD7" / "png" / "PAR_IFN_CD7_A2_P1_layout.png")
    assert img.size == (5 * 10 + 4 * 24, 8)
    assert generate_png_layouts(stats, control, tmp_out, custom) == 0  # kept on re-run
    assert generate_png_layouts(stats, control, tmp_out, custom, overwrite=True) == 3


def test_panels_and_colors():
    assert default_panels(CH) == [("AF488",), ("AF594",), ("DAPI",), ("AF488", "DAPI"), ("AF594", "DAPI")]
    assert channel_color("DAPI") == (0, 0, 1) and channel_color("AF488") == (0, 1, 0)
    assert channel_color("AF594") == (1, 0, 0) and channel_color("Cy5") == (1, 0, 1)
    assert default_panels(["GFP"]) == [("GFP",)]


def test_custom_range_must_be_ordered(czi, tmp_out):
    stats, _ = process_all(czi, PLATE_MAP, SCENES, CH, FIELDS, tmp_out, log=lambda s: None)
    with pytest.raises(ValueError):
        build_control_range(stats, {}, {("CD7", "DAPI"): (10, 10)})


def test_layout_scaling_saturates_cleanly():
    from zeizify.layout import scale_to_uint8

    out = scale_to_uint8(np.array([0, 50, 100, 200]), 50, 100)
    assert out.tolist() == [0, 0, 255, 255]
