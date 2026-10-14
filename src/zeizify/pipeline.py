"""
High-level steps shared by the GUI and the CLI:

    1. map_scenes   - scene -> well mapping + QC outputs (human confirms the QC plot)
    2. extract      - per-well/field/channel MIP TIFFs
    3. finalize     - intensity_report.csv + layout PNGs from the chosen contrast

Each run also writes OUT_DIR/run_manifest.json recording inputs, parameters, software
versions and outputs, so an analysis can be reported and reproduced exactly.
"""

from __future__ import annotations

import json
import platform
import threading
from dataclasses import dataclass, field
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path

from . import __version__
from .contrast import ControlChoice, CustomRanges, write_intensity_report
from .czi_metadata import read_bit_ceiling, read_channel_names, read_field_names, read_pixel_size_um
from .extraction import FieldStat, Log, Progress, process_all
from .layout import Panel, default_color_name, default_panels, generate_png_layouts
from .mapping import SceneWellMapping, build_mapping
from .plate import DEFAULT_PLATE, PlateFormat, load_plate_map
from .qc import write_qc_outputs

QC_SUBDIR = "_scene_mapping"


@dataclass
class CziInfo:
    channel_names: list[str]
    field_names: dict[int, str]
    bit_ceiling: "int | None"
    px_size_um: "float | None"


@dataclass
class RunResult:
    stats: list[FieldStat] = field(default_factory=list)
    skipped: list[tuple[str, str]] = field(default_factory=list)
    info: "CziInfo | None" = None


def read_czi_info(xml: str) -> CziInfo:
    return CziInfo(read_channel_names(xml), read_field_names(xml), read_bit_ceiling(xml), read_pixel_size_um(xml))


def map_scenes(czi_path: Path, out_dir: Path, plate: PlateFormat = DEFAULT_PLATE, flip_rows: bool = False,
               flip_cols: bool = False, log: Log = print) -> tuple[SceneWellMapping, Path]:
    """Build the mapping and write its QC CSV/PNG. Returns (mapping, qc_png_path)."""
    from czifile import CziFile

    qc_dir = out_dir / QC_SUBDIR
    mapping = build_mapping(czi_path, plate, flip_rows, flip_cols, debug_dir=qc_dir)
    log(f"Scene -> well source: {mapping.source}"
        + (" (ZEN ArrayName metadata)" if mapping.source == "metadata" else " (stage-position clustering)"))
    log(f"{len(mapping.scenes)} scenes -> {len(mapping.well_labels)} wells")
    for note in mapping.notes:
        log("  " + note)
    with CziFile(str(czi_path)) as czi:
        field_names = read_field_names(czi.metadata())
    _, png = write_qc_outputs(qc_dir, mapping, field_names)
    return mapping, png


def extract(czi_path: Path, plate_map_path: Path, out_dir: Path, mapping: SceneWellMapping, log: Log = print,
            progress: "Progress | None" = None, cancel: "threading.Event | None" = None) -> RunResult:
    from czifile import CziFile

    plate_map = load_plate_map(plate_map_path)
    if not plate_map:
        raise ValueError(f"No labeled wells in {plate_map_path} - label some wells in the plate map first.")

    with CziFile(str(czi_path)) as czi:
        info = read_czi_info(czi.metadata())
        log(f"Channels: {', '.join(info.channel_names) or '(none found)'}")
        log(f"Sensor ceiling: {info.bit_ceiling} - saturated fields will be flagged" if info.bit_ceiling
            else "Sensor bit depth not in metadata - saturation not checked")
        log(f"Pixel size: {info.px_size_um:.4f} µm/px (embedded in TIFFs)" if info.px_size_um
            else "Pixel size not in metadata - TIFFs carry no physical scale")
        stats, skipped = process_all(czi, plate_map, mapping.well_to_scenes(), info.channel_names, info.field_names,
                                     out_dir, info.bit_ceiling, info.px_size_um, log=log, progress=progress, cancel=cancel)

    sat = sum(s.saturated for s in stats)
    log(f"{len(stats)} channel TIFF(s) under {out_dir}" + (f" - {sat} saturated" if sat else ""))
    for name, reason in skipped:
        log(f"  skipped {name}: {reason}")
    return RunResult(stats, skipped, info)


def finalize(result: RunResult, out_dir: Path, control: ControlChoice, custom: "CustomRanges | None" = None,
             overwrite_png: bool = False, log: Log = print, panels: "list[Panel] | None" = None,
             colors: "dict[str, str] | None" = None) -> tuple[Path, int]:
    report = write_intensity_report(result.stats, control, out_dir, custom)
    n = generate_png_layouts(result.stats, control, out_dir, custom, panels=panels, overwrite=overwrite_png,
                             colors=colors)
    log(f"Wrote {report.name} and {n} layout PNG(s). Raw TIFFs were not modified.")
    return report, n


def _versions() -> dict[str, str]:
    out = {"zeizify": __version__, "python": platform.python_version()}
    for pkg in ("czifile", "tifffile", "imagecodecs", "numpy", "scipy", "pillow", "matplotlib"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def write_manifest(out_dir: Path, czi_path: Path, plate_map_path: Path, plate: PlateFormat, mapping: SceneWellMapping,
                   flip_rows: bool, flip_cols: bool, result: RunResult, control: "ControlChoice | None" = None,
                   custom: "CustomRanges | None" = None, panels: "list[Panel] | None" = None,
                   colors: "dict[str, str] | None" = None) -> Path:
    czi_path = Path(czi_path)
    manifest = {
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "software": _versions(),
        "platform": platform.platform(),
        "inputs": {
            "czi": str(czi_path),
            "czi_size_bytes": czi_path.stat().st_size if czi_path.exists() else None,
            "plate_map": str(plate_map_path),
            "plate_format": plate.name,
        },
        "mapping": {
            "source": mapping.source,
            "flip_rows": flip_rows if mapping.source == "geometry" else None,
            "flip_cols": flip_cols if mapping.source == "geometry" else None,
            "n_scenes": len(mapping.scenes),
            "n_wells": len(mapping.well_labels),
            "notes": mapping.notes,
        },
        "czi": None if result.info is None else {
            "channels": result.info.channel_names,
            "bit_ceiling": result.info.bit_ceiling,
            "pixel_size_um": result.info.px_size_um,
        },
        "projection": "max intensity over Z, per channel, per field",
        "stains": sorted({f"{s.antibody}: {s.channel_name}={s.stain}" for s in result.stats if s.stain}),
        "layout": None if result.info is None else {
            "panels": ["+".join(p) for p in (panels or default_panels(result.info.channel_names))],
            "colors": {c: (colors or {}).get(c) or default_color_name(c) for c in result.info.channel_names},
        },
        "contrast": {
            "control_fields": {ab: {"well": w, "field": f} for ab, (w, f) in (control or {}).items()},
            "custom_ranges": [{"antibody": ab, "channel": ch, "min": lo, "max": hi}
                              for (ab, ch), (lo, hi) in (custom or {}).items()],
        },
        "outputs": {"channel_tiffs": len(result.stats), "skipped": [{"item": n, "reason": r} for n, r in result.skipped]},
    }
    path = out_dir / "run_manifest.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path
