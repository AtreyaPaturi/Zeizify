"""
Plate geometry and the plate-map CSV format shared by the GUI and the pipeline.

A plate map CSV has one row per well with the columns

    Well, Row, Column, CellType, Treatment, Antibody_<fluor> ..., Antibody, Label

There is one Antibody_<fluor> column per fluorophore channel (e.g. Antibody_AF488,
Antibody_AF594), so a well can carry several stains for co-localization. Antibody is
the well's staining panel - the per-channel antibodies joined with '+' (e.g. CD7+HLA) -
and Label is the slugged join of cell type, treatment and panel. Older single-antibody
files (just an Antibody column) are still read. Wells with an empty Label are skipped.
"""

from __future__ import annotations

import csv
import io
import re
import string
from dataclasses import dataclass
from pathlib import Path

PLATE_MAP_COLUMNS = ["Well", "Row", "Column", "CellType", "Treatment", "Antibody", "Label"]
STAIN_PREFIX = "Antibody_"


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


def stain_column(fluor: str) -> str:
    return f"{STAIN_PREFIX}{fluor}"


def panel_name(stains: dict[str, str]) -> str:
    """Staining panel of a well: its per-channel antibodies, in channel order, joined by '+'."""
    return "+".join(v for v in stains.values() if v)


def plate_map_fluorophores(fieldnames: "list[str] | None") -> list[str]:
    """Fluorophore names declared by Antibody_<fluor> columns, in column order."""
    return [c[len(STAIN_PREFIX):] for c in (fieldnames or []) if c.startswith(STAIN_PREFIX) and len(c) > len(STAIN_PREFIX)]


def write_plate_map(path: Path, rows: dict[str, dict[str, str]], plate: PlateFormat = DEFAULT_PLATE,
                    fluorophores: "list[str] | tuple[str, ...]" = ()) -> None:
    """rows: {well: {"CellType", "Treatment", "Antibody_<fluor>" ... or a single "Antibody"}}.
    Missing wells are written empty."""
    with open(path, "w", newline="", encoding="utf-8") as f:
        f.write(plate_map_csv_text(rows, plate, fluorophores))


def plate_map_csv_text(rows: dict[str, dict[str, str]], plate: PlateFormat = DEFAULT_PLATE,
                       fluorophores: "list[str] | tuple[str, ...]" = ()) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, quoting=csv.QUOTE_MINIMAL, lineterminator="\n")
    stain_cols = [stain_column(f) for f in fluorophores]
    w.writerow(PLATE_MAP_COLUMNS[:5] + stain_cols + PLATE_MAP_COLUMNS[5:])
    for well in plate.wells():
        r, c = parse_well(well)
        info = rows.get(well, {})
        ct, tr = info.get("CellType", ""), info.get("Treatment", "")
        stains = [info.get(col, "") for col in stain_cols]
        ab = "+".join(s for s in stains if s) if fluorophores else info.get("Antibody", "")
        w.writerow([well, r, c, ct, tr, *stains, ab, make_label(ct, tr, ab)])
    return buf.getvalue()


def load_plate_map(path: Path) -> dict[str, dict]:
    """Labeled wells only: {well: row}. Each row also gets "Stains": {fluor: antibody} from the
    Antibody_<fluor> columns. Well IDs are normalized (A01 -> A1); every Label is re-slugged so a
    hand-edited CSV can never inject path separators into output names, and is computed from
    the other columns when the file has no Label column."""
    with open(path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        if "Well" not in (reader.fieldnames or []):
            raise ValueError(f"{path} is missing the required 'Well' column")
        fluors = plate_map_fluorophores(reader.fieldnames)
        has_label = "Label" in (reader.fieldnames or [])
        out: dict[str, dict] = {}
        for row in reader:
            row = {k: (v or "").strip() for k, v in row.items() if k is not None}
            if not row.get("Well"):
                continue
            stains = {fl: row.get(stain_column(fl), "") for fl in fluors}
            stains = {fl: ab for fl, ab in stains.items() if ab}
            if not row.get("Antibody") and stains:
                row["Antibody"] = panel_name(stains)
            if not has_label:
                row["Label"] = make_label(row.get("CellType", ""), row.get("Treatment", ""), row.get("Antibody", ""))
            if not row.get("Label"):
                continue
            well = normalize_well(row["Well"])
            row["Well"] = well
            row["Label"] = slug(row["Label"])
            row["Stains"] = stains
            row.setdefault("Antibody", "")
            row.setdefault("Treatment", "")
            row.setdefault("CellType", "")
            if well in out:
                raise ValueError(f"{path}: well {well} appears more than once")
            out[well] = row
    return out
