"""
Persistent GUI state, stored under ~/.zeizify/ (override with ZEIZIFY_HOME):

    plate_state.json  - the plate map being edited (labels + well assignments per category)
    settings.json     - appearance mode and recent runs

Categories are dynamic: cell type, treatment, and one "stain" category per fluorophore the
user declares (e.g. AF488, AF594), so a well can carry several antibodies for co-localization.
"""

from __future__ import annotations

import csv
import json
import os
import string
import uuid
from datetime import datetime
from pathlib import Path

from ..plate import (
    DEFAULT_PLATE,
    PLATE_FORMATS,
    PlateFormat,
    normalize_well,
    plate_map_csv_text,
    plate_map_fluorophores,
    stain_column,
)
from .theme import CATEGORY_PALETTES

FIXED = ["cellType", "treatment"]
FIXED_TITLES = {"cellType": "Cell type", "treatment": "Treatment"}
FIXED_COLUMNS = {"cellType": "CellType", "treatment": "Treatment"}
STAIN = "stain:"
DEFAULT_FLUOROPHORES = ["AF488", "AF594"]
LEGACY_ANTIBODY_FLUOR = "AF594"  # v1.0 single-antibody maps: the antibody-of-interest channel
FLUOR_SUGGESTIONS = ["AF405", "AF488", "AF546", "AF555", "AF568", "AF594", "AF647", "Cy3", "Cy5", "GFP", "mCherry"]


def stain_key(fluor: str) -> str:
    return f"{STAIN}{fluor}"


def is_stain(cat: str) -> bool:
    return cat.startswith(STAIN)


def fluor_of(cat: str) -> str:
    return cat[len(STAIN):]


def home() -> Path:
    p = Path(os.environ.get("ZEIZIFY_HOME", Path.home() / ".zeizify"))
    p.mkdir(parents=True, exist_ok=True)
    return p


def _read_json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001 - missing/corrupt -> start fresh
        return {}


def _write_json(path: Path, data: dict) -> None:
    try:
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps(data, indent=2), encoding="utf-8")
        tmp.replace(path)
    except Exception:  # noqa: BLE001 - persistence is best-effort, never block the UI
        pass


class PlateState:
    """cats[key] = {"labels": [{"id","name","color"}], "byWell": {well: label_id}} for every key in
    categories (cellType, treatment, stain:<fluor> ...)."""

    def __init__(self, path: "Path | None" = None):
        self.path = path or home() / "plate_state.json"
        data = _read_json(self.path)
        self.plate: PlateFormat = PLATE_FORMATS.get(int(data.get("plateWells", 96)), DEFAULT_PLATE)
        self.fluorophores: list[str] = []
        self.cats: dict[str, dict] = {c: {"labels": [], "byWell": {}} for c in FIXED}

        if data.get("version", 1) >= 2:
            self.fluorophores = list(data.get("fluorophores", []))
            for key, cat in data.get("categories", {}).items():
                self.cats[key] = {"labels": list(cat.get("labels", [])), "byWell": dict(cat.get("byWell", {}))}
            for f in self.fluorophores:
                self.cats.setdefault(stain_key(f), {"labels": [], "byWell": {}})
        elif data:  # v1 (Zeizify 1.0/1.1) or the original plate_mapper.py: one "antibody" category
            original = "plateWells" not in data
            for c in FIXED:
                if c in data:
                    self.cats[c] = {"labels": list(data[c].get("labels", [])), "byWell": dict(data[c].get("byWell", {}))}
            self.fluorophores = [LEGACY_ANTIBODY_FLUOR]
            ab = data.get("antibody", {})
            self.cats[stain_key(LEGACY_ANTIBODY_FLUOR)] = {"labels": list(ab.get("labels", [])),
                                                         "byWell": dict(ab.get("byWell", {}))}
            if original:  # the original generated hues were nearly indistinguishable - reassign distinct ones
                for c in self.categories:
                    pal = self.palette(c)
                    for i, lbl in enumerate(self.cats[c]["labels"]):
                        lbl["color"] = pal[i % len(pal)]
        else:
            self.fluorophores = list(DEFAULT_FLUOROPHORES)
            for f in self.fluorophores:
                self.cats[stain_key(f)] = {"labels": [], "byWell": {}}

    def save(self) -> None:
        _write_json(self.path, {"version": 2, "plateWells": self.plate.n_wells, "fluorophores": self.fluorophores,
                                "categories": {c: self.cats[c] for c in self.categories}})

    # categories
    @property
    def categories(self) -> list[str]:
        return FIXED + [stain_key(f) for f in self.fluorophores]

    def title(self, cat: str) -> str:
        return f"Antibody · {fluor_of(cat)}" if is_stain(cat) else FIXED_TITLES[cat]

    def short_title(self, cat: str) -> str:
        return fluor_of(cat) if is_stain(cat) else FIXED_TITLES[cat]

    def badge(self, cat: str) -> str:
        return string.ascii_uppercase[self.categories.index(cat) % 26]

    def palette(self, cat: str) -> list[str]:
        if not is_stain(cat):
            return CATEGORY_PALETTES[cat]
        pal = CATEGORY_PALETTES["antibody"]
        k = (3 * self.fluorophores.index(fluor_of(cat))) % len(pal) if fluor_of(cat) in self.fluorophores else 0
        return pal[k:] + pal[:k]  # each channel starts on a different hue

    def add_fluorophore(self, name: str) -> bool:
        name = name.strip()
        if not name or name.lower() in (f.lower() for f in self.fluorophores):
            return False
        self.fluorophores.append(name)
        self.cats.setdefault(stain_key(name), {"labels": [], "byWell": {}})
        return True

    def remove_fluorophore(self, name: str) -> None:
        self.fluorophores = [f for f in self.fluorophores if f != name]
        self.cats.pop(stain_key(name), None)

    def rename_fluorophore(self, old: str, new: str) -> bool:
        new = new.strip()
        if not new or new == old or new.lower() in (f.lower() for f in self.fluorophores if f != old):
            return False
        self.fluorophores = [new if f == old else f for f in self.fluorophores]
        self.cats[stain_key(new)] = self.cats.pop(stain_key(old), {"labels": [], "byWell": {}})
        return True

    # labels
    def label(self, cat: str, label_id: "str | None") -> "dict | None":
        return next((l for l in self.cats[cat]["labels"] if l["id"] == label_id), None)

    def label_name(self, cat: str, label_id: "str | None") -> str:
        lbl = self.label(cat, label_id)
        return lbl["name"] if lbl else ""

    def name_at(self, cat: str, well: str) -> str:
        return self.label_name(cat, self.cats[cat]["byWell"].get(well))

    def color_at(self, cat: str, well: str) -> "str | None":
        lbl = self.label(cat, self.cats[cat]["byWell"].get(well))
        return lbl["color"] if lbl else None

    def add_label(self, cat: str, name: str) -> dict:
        existing = next((l for l in self.cats[cat]["labels"] if l["name"] == name), None)
        if existing:
            return existing
        used = {l["color"] for l in self.cats[cat]["labels"]}
        pal = self.palette(cat)
        color = next((c for c in pal if c not in used), pal[len(self.cats[cat]["labels"]) % len(pal)])
        lbl = {"id": str(uuid.uuid4()), "name": name, "color": color}
        self.cats[cat]["labels"].append(lbl)
        return lbl

    def delete_label(self, cat: str, label_id: str) -> list[str]:
        self.cats[cat]["labels"] = [l for l in self.cats[cat]["labels"] if l["id"] != label_id]
        changed = [w for w, v in self.cats[cat]["byWell"].items() if v == label_id]
        for w in changed:
            del self.cats[cat]["byWell"][w]
        return changed

    def count(self, cat: str, label_id: str) -> int:
        return sum(1 for v in self.cats[cat]["byWell"].values() if v == label_id)

    # plate
    def set_plate(self, plate: PlateFormat) -> None:
        self.plate = plate
        valid = set(plate.wells())
        for c in self.categories:
            self.cats[c]["byWell"] = {w: v for w, v in self.cats[c]["byWell"].items() if w in valid}

    def clear_wells(self) -> None:
        for c in self.categories:
            self.cats[c]["byWell"] = {}

    def n_labeled(self) -> int:
        return len({w for c in self.categories for w in self.cats[c]["byWell"]})

    # CSV
    def column(self, cat: str) -> str:
        return stain_column(fluor_of(cat)) if is_stain(cat) else FIXED_COLUMNS[cat]

    def rows(self) -> dict[str, dict[str, str]]:
        return {w: {self.column(c): self.name_at(c, w) for c in self.categories} for w in self.plate.wells()}

    def csv_text(self) -> str:
        return plate_map_csv_text(self.rows(), self.plate, self.fluorophores)

    def import_csv(self, path: Path) -> int:
        """Replace the map from a plate map CSV. Fluorophores come from its Antibody_<fluor>
        columns; a single-antibody CSV maps its Antibody column to one fluorophore."""
        with open(path, newline="", encoding="utf-8-sig") as f:
            reader = csv.DictReader(f)
            fields = reader.fieldnames or []
            rows = list(reader)
        wells = [normalize_well(r["Well"]) for r in rows if r.get("Well")]
        fitting = [p for p in PLATE_FORMATS.values() if set(wells) <= set(p.wells())]
        self.plate = min(fitting, key=lambda p: p.n_wells) if fitting else self.plate

        fluors = plate_map_fluorophores(fields)
        columns = {c: FIXED_COLUMNS[c] for c in FIXED}
        if fluors:
            self.fluorophores = fluors
            columns.update({stain_key(f): stain_column(f) for f in fluors})
        else:
            if not self.fluorophores:
                self.fluorophores = [LEGACY_ANTIBODY_FLUOR]
            columns[stain_key(self.fluorophores[0])] = "Antibody"
        self.cats = {c: {"labels": [], "byWell": {}} for c in self.categories}
        for r in rows:
            if not r.get("Well"):
                continue
            w = normalize_well(r["Well"])
            for cat, col in columns.items():
                v = (r.get(col) or "").strip()
                if v:
                    self.cats[cat]["byWell"][w] = self.add_label(cat, v)["id"]
        return self.n_labeled()


class Settings:
    def __init__(self, path: "Path | None" = None):
        self.path = path or home() / "settings.json"
        self.data = _read_json(self.path)
        self.data.setdefault("appearance", "system")
        self.data.setdefault("recents", [])

    def save(self) -> None:
        _write_json(self.path, self.data)

    @property
    def appearance(self) -> str:
        return self.data["appearance"]

    @appearance.setter
    def appearance(self, mode: str) -> None:
        self.data["appearance"] = mode
        self.save()

    @property
    def recents(self) -> list[dict]:
        return self.data["recents"]

    def add_recent(self, czi: str, out: str, plate_map: str) -> None:
        r = [x for x in self.recents if x.get("out") != out]
        r.insert(0, {"czi": czi, "out": out, "plate_map": plate_map, "when": datetime.now().isoformat(timespec="minutes")})
        self.data["recents"] = r[:8]
        self.save()
