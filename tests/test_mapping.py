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
    assert notes == []


@pytest.mark.parametrize("rows,cols", [(True, False), (False, True), (True, True)])
def test_geometry_flips(rows, cols):
    pos, truth = grid_positions()
    got, _ = geometric_mapping(pos, flip_rows=rows, flip_cols=cols)
    assert got == {s: _flip(w, 8, 12, rows, cols) for s, w in truth.items()}


def test_geometry_tolerates_uneven_field_counts():
    pos, truth = grid_positions(skip={"C7", "H12"})
    got, notes = geometric_mapping(pos)
    assert got == truth
    assert "C7=1" in notes[0] and "H12=1" in notes[0]


def test_geometry_384_plate():
    pos, truth = grid_positions(n_rows=16, n_cols=24, fields=2, pitch=4500, field_step=600, tilt=0.002)
    got, _ = geometric_mapping(pos, PLATE_FORMATS[384])
    assert got == truth


def test_geometry_refuses_too_few_scenes():
    pos, _ = grid_positions(n_rows=4, n_cols=6)
    with pytest.raises(MappingError):
        geometric_mapping(pos)


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
    assert any("agrees" in n for n in m.notes)
    assert len(m.well_labels) == 96 and m.well_centroids.shape == (96, 2)


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
