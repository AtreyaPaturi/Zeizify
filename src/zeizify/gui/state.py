"""
Persistent GUI state, stored under ~/.zeizify/ (override with ZEIZIFY_HOME):

    plate_state.json  - the plate map being edited (labels + well assignments per category)
    settings.json     - appearance mode and recent runs
"""

from __future__ import annotations

import csv
import json
import os
import uuid
from datetime import datetime
from pathlib import Path

from ..plate import CATEGORY_COLUMNS, DEFAULT_PLATE, PLATE_FORMATS, PlateFormat, normalize_well, plate_map_csv_text
from .theme import CATEGORY_PALETTES

CATEGORIES = ["cellType", "treatment", "antibody"]
CATEGORY_TITLES = {"cellType": "Cell type", "treatment": "Treatment", "antibody": "Antibody"}
CATEGORY_BADGES = {"cellType": "A", "treatment": "B", "antibody": "C"}


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
    """state[cat] = {"labels": [{"id","name","color"}], "byWell": {well: label_id}}"""

    def __init__(self, path: "Path | None" = None):
        self.path = path or home() / "plate_state.json"
        self.plate: PlateFormat = DEFAULT_PLATE
        self.cats = {c: {"labels": [], "byWell": {}} for c in CATEGORIES}
        data = _read_json(self.path)
        self.plate = PLATE_FORMATS.get(int(data.get("plateWells", 96)), DEFAULT_PLATE)
        legacy = bool(data) and "plateWells" not in data  # written by the original plate_mapper.py
        for c in CATEGORIES:
            if c in data:
                self.cats[c]["labels"] = list(data[c].get("labels", []))
                self.cats[c]["byWell"] = dict(data[c].get("byWell", {}))
            if legacy:  # the old generated hues were nearly indistinguishable - reassign distinct ones
                pal = CATEGORY_PALETTES[c]
                for i, lbl in enumerate(self.cats[c]["labels"]):
                    lbl["color"] = pal[i % len(pal)]

    def save(self) -> None:
        _write_json(self.path, {"plateWells": self.plate.n_wells, **self.cats})

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
        pal = CATEGORY_PALETTES[cat]
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
        for c in CATEGORIES:
            self.cats[c]["byWell"] = {w: v for w, v in self.cats[c]["byWell"].items() if w in valid}

    def clear_wells(self) -> None:
        for c in CATEGORIES:
            self.cats[c]["byWell"] = {}

    def n_labeled(self) -> int:
        return len({w for c in CATEGORIES for w in self.cats[c]["byWell"]})

    # CSV
    def rows(self) -> dict[str, dict[str, str]]:
        out = {}
        for w in self.plate.wells():
            out[w] = {col: self.name_at(cat, w) for cat, col in CATEGORY_COLUMNS.items()}
        return out

    def csv_text(self) -> str:
        return plate_map_csv_text(self.rows(), self.plate)

    def import_csv(self, path: Path) -> int:
        """Replace well assignments from a plate map CSV; labels are created as needed."""
        with open(path, newline="", encoding="utf-8-sig") as f:
            rows = list(csv.DictReader(f))
        wells = [normalize_well(r["Well"]) for r in rows if r.get("Well")]
        n = max((PLATE_FORMATS[k] for k in PLATE_FORMATS if set(wells) <= set(PLATE_FORMATS[k].wells())),
                key=lambda p: -p.n_wells, default=self.plate)
        self.plate = n
        self.clear_wells()
        for r in rows:
            if not r.get("Well"):
                continue
            w = normalize_well(r["Well"])
            for cat, col in CATEGORY_COLUMNS.items():
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
