"""
Scene -> well mapping for a multi-well-plate CZI acquired as scenes = wells x fields.

Primary method - metadata. Every <Scene> element written by ZEN for a plate
acquisition carries an <ArrayName> (e.g. "A1"): the well ID assigned by the
acquisition software itself. This is ground truth, has no orientation ambiguity and
is used whenever every scene has one. When the scene's <Shape> also carries
RowIndex/ColumnIndex, those are cross-checked against ArrayName as an independent
consistency test (the same fields ZEISS' czitools exposes as well_rowID/well_colID).

Fallback method - geometry. If any scene lacks ArrayName (older acquisitions, other
software), well identity is recovered from each scene's absolute stage position (the
subblock start offset on the Y/X axes; metadata only, no pixel data is decoded):
  1. Take one subblock's (Y, X) start per scene.
  2. Agglomeratively cluster scenes into n_wells groups - fields of one well sit far
     closer together than two different wells do. Field counts may differ per well.
  3. Cluster the well centroids' Y values into n_rows bands and, within each band,
     their X values into n_cols bands; every band must come out full or we refuse.
  4. Map band index -> row letter / column number.
Step 4 is a genuine ambiguity (geometry tells us adjacency, not which physical corner
is A1); it is resolved by flip_rows/flip_cols and confirmed by a human on the QC plot.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage

from .czi_metadata import SceneInfo, read_scenes
from .plate import DEFAULT_PLATE, PlateFormat, normalize_well

Log = Callable[[str], None]


class MappingError(ValueError):
    """Scene -> well mapping could not be established with confidence."""

    def __init__(self, message: str, debug_plot: "Path | None" = None):
        super().__init__(message)
        self.debug_plot = debug_plot


@dataclass
class SceneWellMapping:
    scene_to_well: dict[int, str]
    field_points: np.ndarray  # (n_scenes, 2) stage (Y, X), ordered like `scenes`
    scenes: list[int]
    well_labels: np.ndarray  # one per well present
    well_centroids: np.ndarray  # (n_wells, 2) (Y, X)
    source: str  # "metadata" or "geometry"
    notes: list[str] = field(default_factory=list)

    def well_to_scenes(self) -> dict[str, list[int]]:
        out: dict[str, list[int]] = {}
        for scene in sorted(self.scene_to_well):
            out.setdefault(self.scene_to_well[scene], []).append(scene)
        return out


def scene_positions(entries: Iterable) -> dict[int, tuple[int, int]]:
    """One representative (Y, X) stage start per scene index, from subblock directory entries."""
    positions: dict[int, tuple[int, int]] = {}
    for entry in entries:
        s = entry.scene_index
        if s not in positions:
            dims = entry.dims
            positions[s] = (entry.start[dims.index("Y")], entry.start[dims.index("X")])
    return positions


def mapping_from_metadata(scenes: list[SceneInfo]) -> "dict[int, str] | None":
    """Scene -> well from ArrayName. None (never a partial result) if any scene lacks it, so
    the caller falls back to geometry rather than silently trusting an incomplete mapping."""
    if not scenes:
        return None
    mapping: dict[int, str] = {}
    for s in scenes:
        if not s.array_name:
            return None
        try:
            mapping[s.index] = normalize_well(s.array_name)
        except ValueError:
            return None
    return mapping


def crosscheck_shape_indices(scenes: list[SceneInfo], plate: PlateFormat) -> "str | None":
    """Independent check of ArrayName against Shape/RowIndex+ColumnIndex. Returns a human-
    readable note, or None if the indices aren't recorded. Accepts either a 0- or 1-based
    convention, as long as one convention fits every scene."""
    if not scenes or any(s.row_index is None or s.col_index is None for s in scenes):
        return None
    letters = plate.row_letters
    for base in (0, 1):
        ok = True
        for s in scenes:
            r = s.row_index - base
            if not (0 <= r < len(letters)) or normalize_well(s.array_name) != f"{letters[r]}{s.col_index - base + 1}":
                ok = False
                break
        if ok:
            return f"ArrayName agrees with Shape Row/ColumnIndex ({base}-based) for all {len(scenes)} scenes."
    return "WARNING: ArrayName disagrees with Shape Row/ColumnIndex for some scenes - inspect the QC plot."


def cluster_1d(values: np.ndarray, n_clusters: int) -> np.ndarray:
    """Cluster a 1D array into n_clusters groups; labels 0..n-1 ordered by increasing value."""
    if n_clusters == 1:
        return np.zeros(len(values), dtype=int)
    Z = linkage(values.reshape(-1, 1), method="average")
    raw = fcluster(Z, t=n_clusters, criterion="maxclust")  # 1-indexed, arbitrary order
    uniq = np.unique(raw)
    order = np.argsort([values[raw == c].mean() for c in uniq])
    remap = {old: new for new, old in enumerate(uniq[order])}
    return np.array([remap[c] for c in raw])


def geometric_mapping(
    positions: dict[int, tuple[int, int]],
    plate: PlateFormat = DEFAULT_PLATE,
    flip_rows: bool = False,
    flip_cols: bool = False,
    debug_dir: "Path | None" = None,
) -> tuple[dict[int, str], list[str]]:
    """Fallback path: scene -> well purely from stage position. Returns (mapping, notes)."""
    from .qc import save_debug_plot

    scenes = sorted(positions)
    pts = np.array([positions[s] for s in scenes], dtype=float)
    n_rows, n_cols, n_wells = plate.n_rows, plate.n_cols, plate.n_wells
    if len(scenes) < n_wells:
        raise MappingError(
            f"only {len(scenes)} scenes found, fewer than the {n_wells} wells of a {plate.name} plate - "
            f"is the plate format right?"
        )

    # Step 1: fields -> wells.
    field_cluster = fcluster(linkage(pts, method="average"), t=n_wells, criterion="maxclust")
    well_ids = np.unique(field_cluster)
    if len(well_ids) != n_wells:
        sizes = sorted(int(s) for s in np.bincount(field_cluster)[1:] if s > 0)
        raise MappingError(f"expected {n_wells} well clusters, got {len(well_ids)} (sizes {sizes})")
    centroids = np.array([pts[field_cluster == w].mean(axis=0) for w in well_ids])
    field_counts = np.array([int((field_cluster == w).sum()) for w in well_ids])

    def fail(msg: str, row_band: np.ndarray) -> MappingError:
        plot = save_debug_plot(debug_dir, pts, centroids, row_band, msg) if debug_dir else None
        return MappingError(msg + (f" See {plot}." if plot else ""), plot)

    # Step 2/3: snap well centroids onto the grid.
    row_band = cluster_1d(centroids[:, 0], n_rows)
    row_sizes = [int((row_band == r).sum()) for r in range(n_rows)]
    if row_sizes != [n_cols] * n_rows:
        raise fail(f"row-band sizes {row_sizes}, expected {n_cols} in each of {n_rows} rows "
                   f"(a well likely sits near a row boundary - plate tilt?)", row_band)
    col_band = np.empty(n_wells, dtype=int)
    for r in range(n_rows):
        idx = np.where(row_band == r)[0]
        cb = cluster_1d(centroids[idx, 1], n_cols)
        col_sizes = [int((cb == c).sum()) for c in range(n_cols)]
        if col_sizes != [1] * n_cols:
            raise fail(f"row-band {r} column sizes {col_sizes}, expected one well per column", row_band)
        col_band[idx] = cb

    # Step 4: band -> label, with the human-confirmed orientation.
    row_i = (n_rows - 1 - row_band) if flip_rows else row_band
    col_i = (n_cols - 1 - col_band) if flip_cols else col_band
    letters = plate.row_letters
    well_label = np.array([f"{letters[r]}{c + 1}" for r, c in zip(row_i, col_i)])

    notes = []
    modal = int(np.bincount(field_counts).argmax())
    uneven = sorted((str(w), int(c)) for w, c in zip(well_label, field_counts) if c != modal)
    if uneven:
        notes.append(f"{len(uneven)} well(s) differ from the usual {modal} fields: "
                     + ", ".join(f"{w}={c}" for w, c in uneven))

    label_of_cluster = dict(zip(well_ids.tolist(), well_label.tolist()))
    return {s: label_of_cluster[int(c)] for s, c in zip(scenes, field_cluster)}, notes


def _assemble(scene_to_well: dict[int, str], positions: dict[int, tuple[int, int]], source: str, notes: list[str]) -> SceneWellMapping:
    scenes = sorted(positions)
    pts = np.array([positions[s] for s in scenes], dtype=float).reshape(-1, 2)
    labels = sorted(set(scene_to_well.values()), key=lambda w: (w[0], int(w[1:])) if w[1:].isdigit() else (w, 0))
    centroids = np.array([
        pts[[i for i, s in enumerate(scenes) if scene_to_well.get(s) == w]].mean(axis=0) for w in labels
    ]).reshape(-1, 2)
    return SceneWellMapping(scene_to_well, pts, scenes, np.array(labels), centroids, source, notes)


def mapping_from_inputs(
    positions: dict[int, tuple[int, int]],
    scenes_meta: list[SceneInfo],
    plate: PlateFormat = DEFAULT_PLATE,
    flip_rows: bool = False,
    flip_cols: bool = False,
    debug_dir: "Path | None" = None,
) -> SceneWellMapping:
    """Core decision logic, independent of file I/O (see build_mapping)."""
    meta = mapping_from_metadata(scenes_meta)
    if meta is not None and set(positions) <= set(meta):
        scene_to_well = {s: meta[s] for s in positions}
        expected = set(plate.wells())
        unexpected = sorted(set(scene_to_well.values()) - expected)
        if unexpected:
            raise MappingError(
                f"ArrayName references well(s) outside a {plate.name} plate: {unexpected[:10]} - "
                f"choose the right plate format."
            )
        notes = []
        check = crosscheck_shape_indices([s for s in scenes_meta if s.index in positions], plate)
        if check:
            notes.append(check)
        missing = sorted(expected - set(scene_to_well.values()))
        if missing:
            notes.append(f"{len(missing)} well(s) have no scenes: {', '.join(missing[:24])}"
                         + (" ..." if len(missing) > 24 else ""))
        return _assemble(scene_to_well, positions, "metadata", notes)

    scene_to_well, notes = geometric_mapping(positions, plate, flip_rows, flip_cols, debug_dir)
    return _assemble(scene_to_well, positions, "geometry", notes)


def build_mapping(
    czi_path: "str | Path",
    plate: PlateFormat = DEFAULT_PLATE,
    flip_rows: bool = False,
    flip_cols: bool = False,
    debug_dir: "Path | None" = None,
) -> SceneWellMapping:
    """Open a CZI and map its scenes to wells. flip_rows/flip_cols only affect the
    geometry fallback. If debug_dir is given, a failed fallback leaves a diagnostic plot."""
    from czifile import CziFile

    with CziFile(str(czi_path)) as czi:
        positions = scene_positions(czi.filtered_subblock_directory)
        try:
            scenes_meta = read_scenes(czi.metadata())
        except Exception:  # noqa: BLE001 - unreadable metadata just means "use geometry"
            scenes_meta = []
    return mapping_from_inputs(positions, scenes_meta, plate, flip_rows, flip_cols, debug_dir)
