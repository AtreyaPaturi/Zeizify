"""Synthetic fields with known nuclei, spots and colocalization, for the analysis tests."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import tifffile

PX = 0.25  # µm / px -> a 14 µm nucleus is 56 px across
SHAPE = (640, 640)


def _blob(img, y, x, amp, sigma):
    r = int(4 * sigma)
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    g = amp * np.exp(-(yy ** 2 + xx ** 2) / (2 * sigma ** 2))
    img[y - r:y + r + 1, x - r:x + r + 1] += g


def make_field(seed=0, n_spots=(3, 4, 5, 6, 7, 8, 9, 10, 11), shared=0.5, b_extra=2, dapi_scale=None):
    """Nine nuclei on a 3x3 grid. Channel A has n_spots[i] spots in nucleus i; channel B has the first
    round(shared * n) of them at the same place plus b_extra independent spots. Returns
    (dapi, a, b, truth) where truth = list of (n_a, n_shared, n_b) per nucleus in grid order."""
    rng = np.random.default_rng(seed)
    dapi = np.full(SHAPE, 100.0)
    a = np.full(SHAPE, 100.0)
    b = np.full(SHAPE, 100.0)
    yy, xx = np.mgrid[:SHAPE[0], :SHAPE[1]]
    truth = []
    centers = [(110 + 210 * i, 110 + 210 * j) for i in range(3) for j in range(3)]
    for k, (cy, cx) in enumerate(centers):
        ry, rx = 30, 24
        inside = ((yy - cy) / ry) ** 2 + ((xx - cx) / rx) ** 2 <= 1
        dapi[inside] += 2000 * (dapi_scale[k] if dapi_scale else 1)
        a[inside] += 300
        b[inside] += 300
        pts = []
        tries = 0
        while len(pts) < n_spots[k] + b_extra:
            tries += 1
            if tries > 20000:
                raise ValueError(f"cannot fit {n_spots[k] + b_extra} separated spots in one synthetic nucleus")
            py, px_ = cy + rng.integers(-ry + 7, ry - 6), cx + rng.integers(-rx + 7, rx - 6)
            if ((py - cy) / (ry - 6)) ** 2 + ((px_ - cx) / (rx - 6)) ** 2 > 1:
                continue
            if all((py - q[0]) ** 2 + (px_ - q[1]) ** 2 >= 81 for q in pts):
                pts.append((py, px_))
        na = n_spots[k]
        ns = int(round(shared * na))
        for py, px_ in pts[:na]:
            _blob(a, py, px_, 3000, 1.6)
        for py, px_ in pts[:ns] + pts[na:na + b_extra]:
            _blob(b, py, px_, 3000, 1.6)
        truth.append((na, ns, ns + b_extra))
    noise = lambda: rng.normal(0, 15, SHAPE)  # noqa: E731
    to16 = lambda im: np.clip(im + noise(), 0, 65535).astype(np.uint16)  # noqa: E731
    return to16(dapi), to16(a), to16(b), truth


def write_zeizify_output(out: Path, design: list[dict]) -> None:
    """A minimal Zeizify output folder: per-channel TIFFs, intensity_report.csv and run_manifest.json.
    design: [{"well", "field", "treatment", "seed", "shared"}]."""
    out.mkdir(parents=True, exist_ok=True)
    folder = out / "X_Y"
    folder.mkdir(exist_ok=True)
    rows = []
    for d in design:
        dapi, a, b, _ = make_field(seed=d["seed"], n_spots=(6,) * 9, shared=d["shared"])
        label = f"WT_{d['treatment']}_X_Y"
        for ch, img, stain in (("AF488", a, "X"), ("AF594", b, "Y"), ("DAPI", dapi, "")):
            p = folder / f"{label}_{d['well']}_{d['field']}_{ch}.tif"
            tifffile.imwrite(p, img, imagej=True, resolution=(1 / PX, 1 / PX), metadata={"unit": "um"})
            rows.append({"Well": d["well"], "Field": d["field"], "CellType": "WT", "Treatment": d["treatment"],
                         "Antibody": "X+Y", "Label": label, "Channel": ch, "Stain": stain, "OwnMin": int(img.min()),
                         "OwnMax": int(img.max()), "Saturated": False, "ControlWell": "", "ControlField": "",
                         "ControlMin": "", "ControlMax": "", "Path": str(p)})
    with open(out / "intensity_report.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    (out / "run_manifest.json").write_text(json.dumps({"czi": {"pixel_size_um": PX, "bit_ceiling": 16383},
                                                       "inputs": {"plate_map": ""}}), encoding="utf-8")
