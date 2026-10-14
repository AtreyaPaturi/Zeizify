"""
Reference contrast ranges per antibody and the sidecar intensity report.

Fields are grouped by antibody panel - the well's Antibody value, which for co-localization
plate maps is the per-channel stains joined by '+' (e.g. CD7+HLA), so only wells stained
identically are compared. For each panel, one field (typically a mock/vehicle-treated well)
is chosen as the reference; its per-channel min/max becomes the display range for every other field of
that antibody, so treated wells are visually comparable to the control instead of each
auto-scaling to its own signal. Any channel may additionally be given a hand-entered
range per antibody that overrides the reference field for that channel only.

Ranges are written to intensity_report.csv. The raw TIFFs are never modified.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from .extraction import FieldStat

CONTROL_KEYWORDS = ("mock", "control", "ctrl", "dmso", "untreated", "vehicle", "pbs")

ControlChoice = dict[str, tuple[str, str]]  # antibody -> (well, field)
CustomRanges = dict[tuple[str, str], tuple[int, int]]  # (antibody, channel) -> (lo, hi)


def fields_by_antibody(stats: list[FieldStat]) -> dict[str, dict[tuple[str, str], dict]]:
    """{antibody: {(well, field): {"well","field","treatment","channels": {name: (min,max)}}}}"""
    out: dict[str, dict[tuple[str, str], dict]] = defaultdict(dict)
    for s in stats:
        e = out[s.antibody].setdefault(
            (s.well, s.field), {"well": s.well, "field": s.field, "treatment": s.treatment, "channels": {}}
        )
        e["channels"][s.channel_name] = (s.min, s.max)
    return dict(out)


def _well_key(well: str):
    return (well[0], int(well[1:])) if well[1:].isdigit() else (well, 0)


def default_control_choice(stats: list[FieldStat]) -> ControlChoice:
    """First field whose treatment looks like a control (keyword match), else the first field."""
    choice: ControlChoice = {}
    for ab, fields in fields_by_antibody(stats).items():
        entries = sorted(fields.values(), key=lambda e: (_well_key(e["well"]), e["field"]))
        pick = next((e for e in entries if any(k in e["treatment"].lower() for k in CONTROL_KEYWORDS)), entries[0])
        choice[ab] = (pick["well"], pick["field"])
    return choice


def build_control_range(stats: list[FieldStat], control: ControlChoice, custom: "CustomRanges | None" = None) -> dict[tuple[str, str], tuple[int, int]]:
    """{(antibody, channel): (lo, hi)} from each antibody's control field, custom ranges winning."""
    rng: dict[tuple[str, str], tuple[int, int]] = {}
    for s in stats:
        if control.get(s.antibody) == (s.well, s.field):
            rng[(s.antibody, s.channel_name)] = (s.min, s.max)
    for key, (lo, hi) in (custom or {}).items():
        if lo >= hi:
            raise ValueError(f"custom range for {key}: min {lo} must be < max {hi}")
        rng[key] = (lo, hi)
    return rng


REPORT_COLUMNS = [
    "Well", "Field", "CellType", "Treatment", "Antibody", "Label", "Channel", "Stain",
    "OwnMin", "OwnMax", "Saturated", "ControlWell", "ControlField", "ControlMin", "ControlMax", "Path",
]


def write_intensity_report(stats: list[FieldStat], control: ControlChoice, out_dir: Path, custom: "CustomRanges | None" = None) -> Path:
    custom = custom or {}
    rng = build_control_range(stats, control, custom)
    path = out_dir / "intensity_report.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(REPORT_COLUMNS)
        for s in stats:
            if (s.antibody, s.channel_name) in custom:
                cwell, cfield = "custom", ""
            else:
                cwell, cfield = control.get(s.antibody, ("", ""))
            cmin, cmax = rng.get((s.antibody, s.channel_name), ("", ""))
            w.writerow([s.well, s.field, s.cell_type, s.treatment, s.antibody, s.label, s.channel_name, s.stain,
                        s.min, s.max, s.saturated, cwell, cfield, cmin, cmax, str(s.path)])
    return path
