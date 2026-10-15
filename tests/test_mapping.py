import numpy as np
import pytest

from conftest import grid_positions
from zeizify.czi_metadata import SceneInfo
from zeizify.mapping import (
    MappingError,
    cluster_1d,
    crosscheck_shape_indices,
    geometric_mapping,
    mapping_from_inputs,
    mapping_from_metadata,
)
from zeizify.plate import PLATE_FORMATS


def _flip(well: str, n_rows: int, n_cols: int, rows: bool, cols: bool) -> str:
    letters = "ABCDEFGHIJKLMNOP"
    r, c = letters.index(well[0]), int(well[1:]) - 1
    r = n_rows - 1 - r if rows else r
    c = n_cols - 1 - c if cols else c
    return f"{letters[r]}{c + 1}"


def test_cluster_1d_orders_by_value():
    v = np.array([10.0, 500, 11, 1000, 499, 1001])
    assert cluster_1d(v, 3).tolist() == [0, 1, 0, 2, 1, 2]


def test_geometry_recovers_tilted_jittered_plate():
    pos, truth = grid_positions()
    got, notes = geometric_mapping(pos)
    assert got == truth
    assert notes[0].startswith("Detected 96 imaged well(s) in 8 row(s) x 12 column(s), A1–H12")


@pytest.mark.parametrize("rows,cols", [(True, False), (False, True), (True, True)])
def test_geometry_flips(rows, cols):
    pos, truth = grid_positions()
    got, _ = geometric_mapping(pos, flip_rows=rows, flip_cols=cols)
    assert got == {s: _flip(w, 8, 12, rows, cols) for s, w in truth.items()}


def test_geometry_tolerates_uneven_field_counts():
    pos, truth = grid_positions(skip={"C7", "H12"})
    got, notes = geometric_mapping(pos)
    assert got == truth
    assert "C7=1" in notes[-1] and "H12=1" in notes[-1]


def test_geometry_384_plate():
    pos, truth = grid_positions(n_rows=16, n_cols=24, fields=2, pitch=4500, field_step=600, tilt=0.002)
    got, _ = geometric_mapping(pos, PLATE_FORMATS[384])
    assert got == truth


def test_geometry_partial_plate_needs_only_the_first_well():
    pos, truth = grid_positions(n_rows=4, n_cols=6)  # a 4 x 6 block imaged on a 96-well plate
    got, _ = geometric_mapping(pos)
    assert got == truth  # A1-D6
    got, notes = geometric_mapping(pos, first_well="C5")
    shift = {w: f"{chr(ord(w[0]) + 2)}{int(w[1:]) + 4}" for w in truth.values()}
    assert got == {s: shift[w] for s, w in truth.items()}
    assert "C5–F10" in notes[0]
    with pytest.raises(MappingError, match="edge"):
        geometric_mapping(pos, first_well="F8")  # 4 rows from F would run past H


def test_geometry_user_plate_columns_2_to_5():
    """Mirrors a real acquisition: rows A-H x columns 2-5, 5 fields per well spread over most of the
    well, stage pitch ~4% below the nominal 9 mm, pixel size 0.08625 um."""
    import numpy as np

    rng = np.random.default_rng(7)
    px = 0.08625
    pitch = 0.96 * 9000 / px
    offsets = [(-14000, -25000), (0, 0), (14000, 25000), (-14000, 25000), (14000, -25000)]
    pos, truth, s = {}, {}, 0
    for r in range(8):
        for c in range(4):
            for dy, dx in offsets:
                pos[s] = (int(r * pitch + dy + rng.normal(0, 3000)), int(c * pitch + dx + rng.normal(0, 3000)))
                truth[s] = f"{'ABCDEFGH'[r]}{c + 2}"
                s += 1
    got, notes = geometric_mapping(pos, px_size_um=px, first_well="A2")
    assert got == truth
    assert notes[0].startswith("Detected 32 imaged well(s) in 8 row(s) x 4 column(s), A2–H5")


def test_geometry_skipped_column_is_detected_with_pixel_size():
    px = 0.65
    pitch_px = 9000 / px
    pos, truth = grid_positions(n_rows=2, n_cols=4, pitch=int(pitch_px), field_step=900, jitter=30, tilt=0)
    keep = {s: p for s, p in pos.items() if truth[s][1:] != "3"}  # column 3 was not imaged
    got, notes = geometric_mapping(keep, px_size_um=px)
    assert got == {s: truth[s] for s in keep}
    assert "1 column(s) skipped" in notes[1]


def test_geometry_first_well_validation():
    pos, _ = grid_positions(n_rows=2, n_cols=2)
    with pytest.raises(MappingError):
        geometric_mapping(pos, first_well="Z9")


def test_geometry_refuses_ambiguous_rows_and_writes_debug_plot(tmp_path):
    pos, _ = grid_positions(tilt=0.12)  # extreme tilt: rows overlap in Y
    with pytest.raises(MappingError) as ei:
        geometric_mapping(pos, debug_dir=tmp_path)
    assert ei.value.debug_plot is not None and ei.value.debug_plot.exists()


def _scenes_from_truth(truth, base=0, drop=None):
    letters = "ABCDEFGH"
    out = []
    for s, w in truth.items():
        name = "" if s == drop else w
        out.append(SceneInfo(s, "P1", name, letters.index(w[0]) + base, int(w[1:]) - 1 + base))
    return out


def test_metadata_wins_over_geometry_even_if_flips_requested():
    pos, truth = grid_positions()
    m = mapping_from_inputs(pos, _scenes_from_truth(truth), flip_rows=True)
    assert m.source == "metadata"
    assert m.scene_to_well == truth
    assert any("agree with Shape" in n for n in m.notes)
    assert "ArrayName" in m.notes[0]


def test_shape_name_is_used_when_arrayname_is_missing():
    pos, truth = grid_positions()
    letters = "ABCDEFGH"
    scenes = [SceneInfo(s, "P1", "", letters.index(w[0]) + 1, int(w[1:]), shape_name=w) for s, w in truth.items()]
    m = mapping_from_inputs(pos, scenes, flip_rows=True)
    assert m.source == "metadata" and m.scene_to_well == truth
    assert "Shape name" in m.notes[0] and "1-based" in m.notes[1]
    # a non-well shape name (e.g. a tile region) is not trusted -> geometry
    assert len(m.well_labels) == 96 and m.well_centroids.shape == (96, 2)
    scenes[0] = SceneInfo(scenes[0].index, "P1", "", None, None, shape_name="TR1")
    assert mapping_from_inputs(pos, scenes).source == "geometry"


def test_one_missing_arrayname_falls_back_to_geometry():
    pos, truth = grid_positions()
    m = mapping_from_inputs(pos, _scenes_from_truth(truth, drop=5))
    assert m.source == "geometry"
    assert m.scene_to_well == truth
    assert mapping_from_metadata(_scenes_from_truth(truth, drop=5)) is None


def test_crosscheck_detects_disagreement_and_accepts_one_based():
    _, truth = grid_positions()
    assert "1-based" in crosscheck_shape_indices(_scenes_from_truth(truth, base=1), PLATE_FORMATS[96])
    bad = _scenes_from_truth(truth)
    bad[0] = SceneInfo(bad[0].index, "P1", "H12", 0, 0)
    assert crosscheck_shape_indices(bad, PLATE_FORMATS[96]).startswith("WARNING")


def test_metadata_outside_plate_is_an_error():
    pos, truth = grid_positions(n_rows=16, n_cols=24, fields=1, pitch=4500)
    letters = "ABCDEFGHIJKLMNOP"
    scenes = [SceneInfo(s, "P1", w, None, None) for s, w in truth.items()]
    with pytest.raises(MappingError):
        mapping_from_inputs(pos, scenes, PLATE_FORMATS[96])
    assert mapping_from_inputs(pos, scenes, PLATE_FORMATS[384]).source == "metadata"
    assert letters  # silence linters
