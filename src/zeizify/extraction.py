"""
Per-well, per-field, per-channel extraction: read only one scene's subblocks at a time,
max-intensity-project the Z-stack separately per channel, and write each channel once:

    OUT_DIR/<Antibody>/<Label>_<Well>_<Field>_<Channel>.tif

with the CZI's physical pixel size embedded (ImageJ hyperstack metadata, µm). Fields
whose files all exist already are not re-decoded (resumable); their TIFFs are only read
back to recover min/max. Extracted TIFFs are never rewritten by any later step.

Max-intensity projection discards Z by design. Scenes that also contain several mosaic
tiles or time points are refused rather than silently merged into one projection.
"""

from __future__ import annotations

import re
import threading
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import tifffile

from .plate import slug

Log = Callable[[str], None]
Progress = Callable[[int, int], None]


class ExtractionError(RuntimeError):
    pass


@dataclass
class FieldStat:
    well: str
    field: str
    antibody: str
    treatment: str
    cell_type: str
    label: str
    channel: int
    channel_name: str
    min: int
    max: int
    saturated: bool
    path: Path
    stain: str = ""  # antibody in this channel for this well (co-localization plate maps), "" if unknown


def match_stain(channel: str, stains: dict[str, str]) -> str:
    """Antibody assigned to a CZI channel in a well's plate-map stains ({fluor: antibody}).
    Matches the fluorophore name exactly (case/punctuation-insensitive) first, then by its
    wavelength number (e.g. plate map 'Alexa 488' <-> CZI channel 'AF488'), if unambiguous."""
    if not stains:
        return ""
    key = slug(channel).lower()
    for fl, ab in stains.items():
        if slug(fl).lower() == key:
            return ab
    nums = set(re.findall(r"\d{3}", channel))
    hits = [ab for fl, ab in stains.items() if nums & set(re.findall(r"\d{3}", fl))]
    return hits[0] if len(hits) == 1 else ""


def _entry_index(entry, dim: str) -> int:
    return entry.start[entry.dims.index(dim)] if dim in entry.dims else 0


def decode_scene_channel_mips(czi, scene_idx: int, retries: int = 3, retry_delay: float = 2.0, log: Log = print) -> dict[int, np.ndarray]:
    """{channel index: 2D max-intensity projection over Z} for one scene.

    Each subblock read is retried to ride out transient errors on network-mounted files."""
    entries = [e for e in czi.filtered_subblock_directory if e.scene_index == scene_idx]
    if not entries:
        raise ExtractionError(f"no subblocks found for scene {scene_idx}")

    tiles = {(_entry_index(e, "Y"), _entry_index(e, "X"), _entry_index(e, "M")) for e in entries}
    times = {_entry_index(e, "T") for e in entries}
    if len(tiles) > 1 or len(times) > 1:
        raise ExtractionError(
            f"scene {scene_idx} has {len(tiles)} tile position(s) and {len(times)} time point(s); "
            f"Zeizify projects Z only and will not merge tiles or time points"
        )

    by_channel: dict[int, list[np.ndarray]] = defaultdict(list)
    for e in entries:
        c = _entry_index(e, "C")
        for attempt in range(1, retries + 1):
            try:
                plane = np.squeeze(np.asarray(e.asimage(czi).asarray()))
                break
            except OSError as exc:
                if attempt == retries:
                    raise
                log(f"    read error, retrying {attempt}/{retries - 1}: {exc}")
                time.sleep(retry_delay)
        if plane.ndim != 2:
            raise ExtractionError(f"scene {scene_idx} channel {c}: expected 2D planes, got shape {plane.shape}")
        by_channel[c].append(plane)

    return {c: np.max(np.stack(planes), axis=0) for c, planes in by_channel.items()}


def channel_name(channel_names: list[str], c: int) -> str:
    return slug(channel_names[c]) if c < len(channel_names) and slug(channel_names[c]) else f"C{c}"


def output_path(out_dir: Path, row: dict, well: str, field: str, cname: str) -> Path:
    return out_dir / (slug(row.get("Antibody", "")) or "unassigned") / f"{row['Label']}_{well}_{field}_{cname}.tif"


def write_channel_tiff(path: Path, plane: np.ndarray, px_size_um: "float | None") -> None:
    kwargs: dict = {"photometric": "minisblack"}
    if px_size_um:
        kwargs.update(imagej=True, resolution=(1.0 / px_size_um, 1.0 / px_size_um), metadata={"unit": "um"})
    tifffile.imwrite(path, plane, **kwargs)


def process_all(
    czi,
    plate_map: dict[str, dict[str, str]],
    scenes_by_well: dict[str, list[int]],
    channel_names: list[str],
    field_names: dict[int, str],
    out_dir: Path,
    bit_ceiling: "int | None" = None,
    px_size_um: "float | None" = None,
    log: Log = print,
    progress: "Progress | None" = None,
    cancel: "threading.Event | None" = None,
) -> tuple[list[FieldStat], list[tuple[str, str]]]:
    """Extract every labeled well. Returns (stats, skipped) - one bad field never aborts the run."""
    stats: list[FieldStat] = []
    skipped: list[tuple[str, str]] = []
    wells = sorted(plate_map, key=lambda w: (w[0], int(w[1:])))
    total_fields = sum(len(scenes_by_well.get(w, [])) for w in wells)
    done = 0

    def stat(row, well, field, c, cname, lo, hi, path) -> FieldStat:
        return FieldStat(well, field, row.get("Antibody", ""), row.get("Treatment", ""), row.get("CellType", ""),
                         row["Label"], c, cname, lo, hi, bool(bit_ceiling) and hi >= bit_ceiling, path,
                         match_stain(cname, row.get("Stains") or {}))

    for i, well in enumerate(wells, start=1):
        row = plate_map[well]
        scenes = scenes_by_well.get(well)
        if not scenes:
            skipped.append((well, "no scenes mapped to this well"))
            continue
        log(f"[{i}/{len(wells)}] {well} ({row['Label']}) - {len(scenes)} field(s)")

        for order, scene in enumerate(scenes, start=1):
            if cancel is not None and cancel.is_set():
                log("Cancelled - finished files are kept; re-running resumes from here.")
                return stats, skipped
            field = slug(field_names.get(scene, "")) or f"F{order}"
            expected = {c: output_path(out_dir, row, well, field, channel_name(channel_names, c))
                        for c in range(len(channel_names))}

            if expected and all(p.exists() for p in expected.values()):
                log(f"    {field}: already extracted - reusing")
                for c, p in expected.items():
                    arr = tifffile.imread(p)
                    stats.append(stat(row, well, field, c, channel_name(channel_names, c), int(arr.min()), int(arr.max()), p))
            else:
                try:
                    mips = decode_scene_channel_mips(czi, scene, log=log)
                except Exception as exc:  # noqa: BLE001 - report at the end, keep going
                    skipped.append((f"{well} {field} (scene {scene})", str(exc)))
                    done += 1
                    if progress:
                        progress(done, total_fields)
                    continue
                for c in sorted(mips):
                    plane = mips[c]
                    cname = channel_name(channel_names, c)
                    path = output_path(out_dir, row, well, field, cname)
                    path.parent.mkdir(parents=True, exist_ok=True)
                    write_channel_tiff(path, plane, px_size_um)
                    lo, hi = int(plane.min()), int(plane.max())
                    s = stat(row, well, field, c, cname, lo, hi, path)
                    log(f"    {field} {cname}: min={lo} max={hi}" + ("  [SATURATED]" if s.saturated else ""))
                    stats.append(s)
            done += 1
            if progress:
                progress(done, total_fields)

    return stats, skipped
