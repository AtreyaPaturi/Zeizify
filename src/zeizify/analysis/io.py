"""Discover the fields Zeizify extracted, with their channels, stains and plate-map metadata."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path

import tifffile

from ..layout import nuclear_channel
from ..plate import load_plate_map, normalize_well


@dataclass
class FieldImages:
    well: str
    field: str
    antibody: str  # staining panel, e.g. "SUMO2/3+PML"
    treatment: str
    cell_type: str
    label: str
    replicate: str
    channels: dict[str, Path]  # channel name -> TIFF
    stains: dict[str, str] = field(default_factory=dict)  # channel name -> antibody in it
    saturated: dict[str, bool] = field(default_factory=dict)

    @property
    def key(self) -> str:
        return f"{self.label}_{self.well}_{self.field}"


@dataclass
class Dataset:
    out_dir: Path
    fields: list[FieldImages]
    channels: list[str]
    px_um: "float | None"
    bit_ceiling: "int | None"

    @property
    def nuclear(self) -> "str | None":
        return nuclear_channel(self.channels)

    @property
    def signal_channels(self) -> list[str]:
        return [c for c in self.channels if c != self.nuclear]

    def stains_for(self, channel: str) -> list[str]:
        return sorted({f.stains.get(channel, "") for f in self.fields} - {""})


def _resolve(path: str, out_dir: Path, antibody_dir: str) -> Path:
    """The report stores absolute paths; fall back to the same file name under out_dir when the
    output folder was moved or opened from another computer."""
    p = Path(path)
    if p.exists():
        return p
    for cand in (out_dir / antibody_dir / p.name, *out_dir.rglob(p.name)):
        if cand.exists():
            return cand
    return p


def _plate_map_extra(out_dir: Path) -> dict[str, dict]:
    """Plate-map rows (for a 'Replicate' column), from the manifest's plate map or out_dir/plate_map.csv."""
    cands = [out_dir / "plate_map.csv"]
    try:
        m = json.loads((out_dir / "run_manifest.json").read_text(encoding="utf-8"))
        cands.insert(0, Path(m["inputs"]["plate_map"]))
    except Exception:  # noqa: BLE001
        pass
    for c in cands:
        try:
            if c.exists():
                return load_plate_map(c)
        except Exception:  # noqa: BLE001
            continue
    return {}


def load_dataset(out_dir: "str | Path", replicate_column: str = "") -> Dataset:
    from ..plate import slug

    out_dir = Path(out_dir)
    report = out_dir / "intensity_report.csv"
    if not report.exists():
        raise FileNotFoundError(f"{report} not found - run Zeizify's extraction and Apply & export first.")

    pm = _plate_map_extra(out_dir)
    by_field = replicate_column.lower() == "field"  # technical replicates: each field of a well
    if replicate_column.lower() in ("well", "field"):
        rep_col = ""
    else:
        rep_col = replicate_column or ("Replicate" if any(r.get("Replicate") for r in pm.values()) else "")

    fields: dict[tuple[str, str], FieldImages] = {}
    channels: list[str] = []
    with open(report, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            well = normalize_well(r["Well"])
            key = (well, r["Field"])
            rep = (pm.get(well, {}).get(rep_col) or "").strip() if rep_col else ""
            fi = fields.setdefault(key, FieldImages(
                well=well, field=r["Field"], antibody=r.get("Antibody", ""), treatment=r.get("Treatment", ""),
                cell_type=r.get("CellType", ""), label=r.get("Label", ""),
                replicate=f"{well}:{r['Field']}" if by_field else (rep or well), channels={}))
            ch = r["Channel"]
            if ch not in channels:
                channels.append(ch)
            fi.channels[ch] = _resolve(r["Path"], out_dir, slug(fi.antibody) or "unassigned")
            if r.get("Stain"):
                fi.stains[ch] = r["Stain"]
            fi.saturated[ch] = str(r.get("Saturated", "")).lower() == "true"

    px = bits = None
    try:
        m = json.loads((out_dir / "run_manifest.json").read_text(encoding="utf-8"))
        px = (m.get("czi") or {}).get("pixel_size_um")
        bits = (m.get("czi") or {}).get("bit_ceiling")
    except Exception:  # noqa: BLE001
        pass
    if not px and fields:  # fall back to the ImageJ resolution tag written by Zeizify
        try:
            first = next(iter(next(iter(fields.values())).channels.values()))
            with tifffile.TiffFile(first) as t:
                num, den = t.pages[0].tags["XResolution"].value
                px = den / num if num else None
        except Exception:  # noqa: BLE001
            px = None

    ordered = sorted(fields.values(), key=lambda fi: (fi.well[0], int(fi.well[1:]), fi.field))
    return Dataset(out_dir, ordered, channels, px, bits)
