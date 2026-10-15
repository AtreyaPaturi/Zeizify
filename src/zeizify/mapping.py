"""
Scene -> well mapping for a multi-well-plate CZI acquired as scenes = wells x fields.

Primary method - metadata. ZEN records each scene's well in its <Scene> element: as
<ArrayName> (e.g. "A1") in plate-carrier acquisitions, or as the name of the scene's
<Shape> (e.g. <Shape Name="A2">) in others - the same two fields, in the same order of
preference, that ZEISS' czitools reads. This is ground truth, has no orientation ambiguity
and is used whenever every scene has a well ID. When the Shape also carries RowIndex /
ColumnIndex, those are cross-checked against the well ID as an independent consistency test.

Fallback method - geometry. If any scene lacks a well ID, wells are recovered from each
scene's stage position (the subblock start offset on the Y/X axes; metadata only, no pixel
data is decoded). Only part of the plate may have been imaged (e.g. columns 2-5), so the
number of wells is not assumed:
  1. Take one subblock's (Y, X) start per scene.
  2. Split the scenes' Y values into row bands and their X values into column bands at the
     gaps between wells: gaps wider than ~1/3 of the well pitch (from the plate format and
     pixel size), or - when the pixel size is unknown - the clearly-largest gaps.
  3. Check every band is narrower than one well pitch, so that two rows (or columns) can
     never be merged silently, e.g. by plate tilt.
  4. Number the bands by their distance from the first band in units of the well pitch, so
     rows/columns skipped during acquisition are accounted for.
  5. Label wells relative to the first imaged well (default A1, e.g. A2 when imaging started
     in column 2), optionally reversing row/column order.
Steps 4-5 are a genuine ambiguity - geometry tells us adjacency, not which physical well
is which - so the first well and flip_rows/flip_cols are confirmed by a human on the QC plot.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

import numpy as np

from .czi_metadata import SceneInfo, read_pixel_size_um, read_scenes
from .plate import DEFAULT_PLATE, PlateFormat, normalize_well, parse_well

Log = Callable[[str], None]

_WELL_SHAPE = re.compile(r"^[A-Pa-p]0?\d{1,2}$")

# Center-to-center well spacing of ANSI/SLAS (SBS) plates, in mm.
SBS_PITCH_MM = {6: 39.12, 12: 26.01, 24: 19.3, 48: 13.08, 96: 9.0, 384: 4.5}


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
    first_well: "str | None" = None  # geometry only: the imaged well closest to A1

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
    """Scene -> well from ArrayName, else Shape/@Name. None (never a partial result) if any scene
    lacks a valid well ID, so the caller falls back to geometry rather than silently trusting an
    incomplete mapping."""
    if not scenes:
        return None
    mapping: dict[int, str] = {}
    for s in scenes:
        if not s.well_name:
            return None
        # Shape names are free text in non-plate experiments (e.g. tile region "TR1"): only accept the
        # well-ID shape A1..P24. ArrayName is written by ZEN's plate carrier and is trusted as-is.
        if not s.array_name and not _WELL_SHAPE.match(s.shape_name):
            return None
        try:
            mapping[s.index] = normalize_well(s.well_name)
        except ValueError:
            return None
    return mapping


def metadata_tag(scenes: list[SceneInfo]) -> str:
    if all(s.array_name for s in scenes):
        return "ArrayName"
    if not any(s.array_name for s in scenes):
        return "Shape name"
    return "ArrayName / Shape name"


def crosscheck_shape_indices(scenes: list[SceneInfo], plate: PlateFormat) -> "str | None":
    """Independent check of the well ID against Shape/RowIndex+ColumnIndex. Returns a human-
    readable note, or None if the indices aren't recorded. Accepts either a 0- or 1-based
    convention, as long as one convention fits every scene."""
    if not scenes or any(s.row_index is None or s.col_index is None for s in scenes):
        return None
    letters = plate.row_letters
    for base in (0, 1):
        ok = True
        for s in scenes:
            r = s.row_index - base
            if not (0 <= r < len(letters)) or normalize_well(s.well_name) != f"{letters[r]}{s.col_index - base + 1}":
                ok = False
                break
        if ok:
            return f"Well IDs agree with Shape Row/ColumnIndex ({base}-based) for all {len(scenes)} scenes."
    return "WARNING: well IDs disagree with Shape Row/ColumnIndex for some scenes - inspect the QC plot."


def cluster_1d(values: np.ndarray, n_clusters: int) -> np.ndarray:
    """Cluster a 1D array into n_clusters groups; labels 0..n-1 ordered by increasing value."""
    from scipy.cluster.hierarchy import fcluster, linkage

    if n_clusters == 1:
        return np.zeros(len(values), dtype=int)
    Z = linkage(values.reshape(-1, 1), method="average")
    raw = fcluster(Z, t=n_clusters, criterion="maxclust")  # 1-indexed, arbitrary order
    uniq = np.unique(raw)
    order = np.argsort([values[raw == c].mean() for c in uniq])
    remap = {old: new for new, old in enumerate(uniq[order])}
    return np.array([remap[c] for c in raw])


def split_bands(values: np.ndarray, gap: "float | None") -> np.ndarray:
    """Band label (0..k-1, by increasing value) per value, splitting sorted values at gaps > gap.
    With gap=None the threshold is half the largest gap, used only if that gap clearly stands
    out (> 4x the median gap); otherwise everything is one band."""
    order = np.argsort(values)
    d = np.diff(values[order])
    if len(d) == 0:
        return np.zeros(len(values), dtype=int)
    if gap is None:
        md, med = float(d.max()), float(np.median(d))
        gap = md / 2 if md > 4 * max(med, 1.0) else np.inf
    labels_sorted = np.concatenate([[0], np.cumsum(d > gap)])
    labels = np.empty(len(values), dtype=int)
    labels[order] = labels_sorted
    return labels


def _midranges(values: np.ndarray, labels: np.ndarray) -> np.ndarray:
    return np.array([(values[labels == b].min() + values[labels == b].max()) / 2 for b in range(labels.max() + 1)])


def _extents(values: np.ndarray, labels: np.ndarray) -> np.ndarray:
    return np.array([np.ptp(values[labels == b]) for b in range(labels.max() + 1)])


def geometric_mapping(
    positions: dict[int, tuple[int, int]],
    plate: PlateFormat = DEFAULT_PLATE,
    flip_rows: bool = False,
    flip_cols: bool = False,
    debug_dir: "Path | None" = None,
    px_size_um: "float | None" = None,
    first_well: str = "A1",
) -> tuple[dict[int, str], list[str]]:
    """Fallback path: scene -> well purely from stage position. Returns (mapping, notes).
    Works for full and partially imaged plates; first_well is the imaged well closest to A1."""
    from .qc import save_debug_plot

    scenes = sorted(positions)
    pts = np.array([positions[s] for s in scenes], dtype=float).reshape(-1, 2)
    if len(scenes) == 0:
        raise MappingError("no scenes found in the CZI")
    ys, xs = pts[:, 0], pts[:, 1]

    nominal = None
    if px_size_um and plate.n_wells in SBS_PITCH_MM:
        nominal = SBS_PITCH_MM[plate.n_wells] * 1000.0 / px_size_um  # well pitch in pixels

    def fail(msg: str, labels: np.ndarray) -> MappingError:
        plot = save_debug_plot(debug_dir, pts, pts, labels, msg) if debug_dir else None
        return MappingError(msg + (f" See {plot}." if plot else ""), plot)

    # Step 2: bands at the gaps between wells.
    gap = 0.35 * nominal if nominal else None
    row_band, col_band = split_bands(ys, gap), split_bands(xs, gap)
    row_c, col_c = _midranges(ys, row_band), _midranges(xs, col_band)

    # Well pitch actually used: the median spacing between neighbouring bands that are one well apart
    # (stage scaling can differ a few % from nominal, and band centres shift with where the fields
    # happen to sit, so neither the nominal value nor the single closest pair is reliable on its own).
    spacings = np.concatenate([np.diff(row_c), np.diff(col_c)])
    if nominal:
        single = spacings[(spacings > 0.6 * nominal) & (spacings < 1.4 * nominal)]
        pitch = float(np.median(single)) if len(single) else nominal
    elif len(spacings):
        single = spacings[spacings < 1.5 * spacings.min()]
        pitch = float(np.median(single))
    else:
        pitch = None

    # Step 3: no band may be wider than a well.
    if pitch:
        for axis, vals, lab in (("row", ys, row_band), ("column", xs, col_band)):
            ext = _extents(vals, lab)
            if ext.max() > 0.85 * pitch:
                raise fail(f"a {axis} band spans {ext.max() / pitch:.1f} well pitches - fields of neighbouring "
                           f"{axis}s overlap (plate tilt, or the wrong plate format?)", row_band)

    # Step 4: band -> grid index, leaving room for skipped rows/columns.
    def indices(centers: np.ndarray, axis: str) -> np.ndarray:
        if not pitch or len(centers) == 1:
            return np.arange(len(centers))
        steps = (centers - centers[0]) / pitch
        idx = np.rint(steps).astype(int)
        worst = float(np.max(np.abs(steps - idx)))
        if worst > 0.35:
            raise fail(f"{axis} spacing is irregular: a {axis} sits {worst:.2f} well pitches off the grid - cannot "
                       f"tell whether a {axis} was skipped (wrong plate format?)", row_band)
        if np.any(np.diff(idx) < 1):
            raise fail(f"two {axis} bands are closer than one well pitch - cannot tell the {axis}s apart", row_band)
        return idx

    r_idx, c_idx = indices(row_c, "row")[row_band], indices(col_c, "column")[col_band]
    n_r, n_c = int(r_idx.max()) + 1, int(c_idx.max()) + 1
    if flip_rows:
        r_idx = n_r - 1 - r_idx
    if flip_cols:
        c_idx = n_c - 1 - c_idx

    # Step 5: absolute labels from the first imaged well.
    try:
        fr, fc = parse_well(first_well)
    except ValueError as exc:
        raise MappingError(f"first well {first_well!r} is not a well ID") from exc
    letters = plate.row_letters
    if fr not in letters or not 1 <= fc <= plate.n_cols:
        raise MappingError(f"first well {first_well} is not on a {plate.name} plate")
    r0, c0 = letters.index(fr), fc - 1
    if r0 + n_r > plate.n_rows or c0 + n_c > plate.n_cols:
        last = f"{letters[min(r0 + n_r, plate.n_rows) - 1]}{min(c0 + n_c, plate.n_cols)}"
        raise fail(f"the imaged area spans {n_r} rows x {n_c} columns; starting at {fr}{fc} it would run past "
                   f"the edge of a {plate.name} plate (last well {last}). Check the first well and plate format.",
                   row_band)
    wells = [f"{letters[r0 + r]}{c0 + c + 1}" for r, c in zip(r_idx, c_idx)]
    scene_to_well = dict(zip(scenes, wells))

    present = sorted(set(wells), key=lambda w: (w[0], int(w[1:])))
    notes = [f"Detected {len(present)} imaged well(s) in {len(row_c)} row(s) x {len(col_c)} column(s), "
             f"{present[0]}–{present[-1]}" + (f" (well pitch {pitch:.0f} px)" if pitch else "") + "."]
    skipped_r = n_r - len(row_c)
    skipped_c = n_c - len(col_c)
    if skipped_r or skipped_c:
        notes.append(f"Gaps in the imaged grid: {skipped_r} row(s) and {skipped_c} column(s) skipped.")
    counts = {w: wells.count(w) for w in present}
    modal = max(set(counts.values()), key=list(counts.values()).count)
    uneven = [(w, c) for w, c in counts.items() if c != modal]
    if uneven:
        notes.append(f"{len(uneven)} well(s) differ from the usual {modal} fields: "
                     + ", ".join(f"{w}={c}" for w, c in uneven))
    return scene_to_well, notes


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
    px_size_um: "float | None" = None,
    first_well: str = "A1",
) -> SceneWellMapping:
    """Core decision logic, independent of file I/O (see build_mapping)."""
    meta = mapping_from_metadata(scenes_meta)
    if meta is not None and set(positions) <= set(meta):
        scene_to_well = {s: meta[s] for s in positions}
        expected = set(plate.wells())
        unexpected = sorted(set(scene_to_well.values()) - expected)
        if unexpected:
            raise MappingError(
                f"Well IDs in the metadata lie outside a {plate.name} plate: {unexpected[:10]} - "
                f"choose the right plate format."
            )
        used = [s for s in scenes_meta if s.index in positions]
        notes = [f"Well IDs read from each scene's {metadata_tag(used)}."]
        check = crosscheck_shape_indices(used, plate)
        if check:
            notes.append(check)
        present = sorted(set(scene_to_well.values()), key=lambda w: (w[0], int(w[1:])))
        notes.append(f"{len(present)} imaged well(s): {', '.join(present[:24])}" + (" ..." if len(present) > 24 else ""))
        return _assemble(scene_to_well, positions, "metadata", notes)

    scene_to_well, notes = geometric_mapping(positions, plate, flip_rows, flip_cols, debug_dir, px_size_um, first_well)
    mapping = _assemble(scene_to_well, positions, "geometry", notes)
    mapping.first_well = normalize_well(first_well)
    return mapping


def build_mapping(
    czi_path: "str | Path",
    plate: PlateFormat = DEFAULT_PLATE,
    flip_rows: bool = False,
    flip_cols: bool = False,
    debug_dir: "Path | None" = None,
    first_well: str = "A1",
) -> SceneWellMapping:
    """Open a CZI and map its scenes to wells. flip_rows/flip_cols/first_well only affect the
    geometry fallback. If debug_dir is given, a failed fallback leaves a diagnostic plot."""
    from czifile import CziFile

    with CziFile(str(czi_path)) as czi:
        positions = scene_positions(czi.filtered_subblock_directory)
        try:
            xml = czi.metadata()
            scenes_meta, px = read_scenes(xml), read_pixel_size_um(xml)
        except Exception:  # noqa: BLE001 - unreadable metadata just means "use geometry"
            scenes_meta, px = [], None
    return mapping_from_inputs(positions, scenes_meta, plate, flip_rows, flip_cols, debug_dir, px, first_well)
