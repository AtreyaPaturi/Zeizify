"""
Plate geometry and the plate-map CSV format shared by the GUI and the pipeline.

A plate map CSV has one row per well with the columns
    Well, Row, Column, CellType, Treatment, Antibody, Label
where Label is the slugged join of the three category values. Wells with an empty
Label are treated as unlabeled and skipped by the extraction step.
"""

from __future__ import annotations

import csv
import re
import string
from dataclasses import dataclass
from pathlib import Path

PLATE_MAP_COLUMNS = ["Well", "Row", "Column", "CellType", "Treatment", "Antibody", "Label"]
CATEGORY_COLUMNS = {"cellType": "CellType", "treatment": "Treatment", "antibody": "Antibody"}


@dataclass(frozen=True)
class PlateFormat:
    """An SBS-style rectangular plate: rows are letters, columns are 1-based numbers."""

    n_rows: int
    n_cols: int

    @property
    def n_wells(self) -> int:
        return self.n_rows * self.n_cols

    @property
    def row_letters(self) -> str:
        return string.ascii_uppercase[: self.n_rows]

    @property
    def name(self) -> str:
        return f"{self.n_wells}-well"

    def wells(self) -> list[str]:
        """All well IDs in row-major order (A1, A2, ..., B1, ...)."""
        return [f"{r}{c}" for r in self.row_letters for c in range(1, self.n_cols + 1)]

    def row_wells(self, row: str) -> list[str]:
        return [f"{row}{c}" for c in range(1, self.n_cols + 1)]

    def col_wells(self, col: int) -> list[str]:
        return [f"{r}{col}" for r in self.row_letters]


PLATE_FORMATS = {
    6: PlateFormat(2, 3),
    12: PlateFormat(3, 4),
    24: PlateFormat(4, 6),
    48: PlateFormat(6, 8),
    96: PlateFormat(8, 12),
    384: PlateFormat(16, 24),
}
DEFAULT_PLATE = PLATE_FORMATS[96]

_WELL_RE = re.compile(r"^([A-Za-z]+)0*(\d+)$")


def parse_well(well: str) -> tuple[str, int]:
    """'B07' or 'b7' -> ('B', 7). Raises ValueError on anything else."""
    m = _WELL_RE.match(well.strip())
    if not m:
        raise ValueError(f"not a well ID: {well!r}")
    return m.group(1).upper(), int(m.group(2))


def normalize_well(well: str) -> str:
    """Canonical well ID without zero padding: 'a01' -> 'A1'."""
    row, col = parse_well(well)
    return f"{row}{col}"


def slug(s: str) -> str:
    """Filesystem-safe token: runs of non-alphanumerics collapse to '_'."""
    return re.sub(r"[^a-zA-Z0-9]+", "_", (s or "").strip()).strip("_")


def make_label(cell_type: str, treatment: str, antibody: str) -> str:
    return "_".join(p for p in (slug(cell_type), slug(treatment), slug(antibody)) if p)


def write_plate_map(path: Path, rows: dict[str, dict[str, str]], plate: PlateFormat = DEFAULT_PLATE) -> None:
    """rows: {well: {"CellType":..., "Treatment":..., "Antibody":...}} - missing wells are written empty."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(plate_map_csv_text(rows, plate))


def plate_map_csv_text(rows: dict[str, dict[str, str]], plate: PlateFormat = DEFAULT_PLATE) -> str:
    import io

    buf = io.StringIO()
    w = csv.writer(buf, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    w.writerow(PLATE_MAP_COLUMNS)
    for well in plate.wells():
        r, c = parse_well(well)
        info = rows.get(well, {})
        ct, tr, ab = info.get("CellType", ""), info.get("Treatment", ""), info.get("Antibody", "")
        w.writerow([well, r, c, ct, tr, ab, make_label(ct, tr, ab)])
    return buf.getvalue()


def load_plate_map(path: Path) -> dict[str, dict[str, str]]:
    """Labeled wells only: {well: row}. Well IDs are normalized (A01 -> A1) and every Label
    is re-slugged so a hand-edited CSV can never inject path separators into output names."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        missing = {"Well", "Label"} - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing required column(s): {sorted(missing)}")
        out: dict[str, dict[str, str]] = {}
        for row in reader:
            if not (row.get("Label") or "").strip():
                continue
            well = normalize_well(row["Well"])
            row = {k: (v or "").strip() for k, v in row.items() if k is not None}
            row["Well"] = well
            row["Label"] = slug(row["Label"])
            row.setdefault("Antibody", "")
            row.setdefault("Treatment", "")
            row.setdefault("CellType", "")
            if well in out:
                raise ValueError(f"{path}: well {well} appears more than once")
            out[well] = row
    return out
