"""Synthetic stand-ins for CZI files so the pipeline can be tested without real data."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest


@dataclass
class FakeEntry:
    scene_index: int
    dims: tuple
    start: tuple
    data: np.ndarray

    def asimage(self, czi):
        entry = self

        class _Img:
            def asarray(self):
                return entry.data

        return _Img()


class FakeCzi:
    def __init__(self, entries):
        self.filtered_subblock_directory = tuple(entries)


def make_entries(scene_wells: dict[int, tuple[int, int]], n_c=3, n_z=4, shape=(8, 10), seed=0):
    """One entry per (scene, channel, z). Pixel value encodes scene/channel/z so the MIP is predictable:
    the max over z is (scene*100 + c*10 + n_z-1) everywhere, plus a single bright pixel."""
    rng = np.random.default_rng(seed)
    entries = []
    for s, (y, x) in scene_wells.items():
        for c in range(n_c):
            for z in range(n_z):
                arr = np.full((1, 1, *shape, 1), s * 100 + c * 10 + z, dtype=np.uint16)
                if z == 1:
                    arr[0, 0, 0, 0, 0] = 5000 + c
                entries.append(FakeEntry(s, ("C", "Z", "Y", "X", "S"), (c, z, y, x, 0), arr))
    rng.shuffle(entries)
    return entries


def grid_positions(n_rows=8, n_cols=12, fields=3, pitch=9000, field_step=800, jitter=40, tilt=0.004, seed=1,
                   skip=()):
    """Stage (Y, X) per scene for a plate scanned row-major, plus ground-truth wells."""
    rng = np.random.default_rng(seed)
    letters = "ABCDEFGHIJKLMNOP"
    positions, truth = {}, {}
    s = 0
    for r in range(n_rows):
        for c in range(n_cols):
            well = f"{letters[r]}{c + 1}"
            nf = 1 if well in skip else fields
            for f in range(nf):
                y = r * pitch + c * pitch * tilt + rng.normal(0, jitter)
                x = c * pitch + f * field_step + rng.normal(0, jitter)
                positions[s] = (int(y), int(x))
                truth[s] = well
                s += 1
    return positions, truth


@pytest.fixture
def tmp_out(tmp_path):
    d = tmp_path / "out"
    d.mkdir()
    return d
