"""Quality-control outputs for the scene -> well mapping (CSV + stage-position plot)."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from .mapping import SceneWellMapping


def _pyplot():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def save_debug_plot(debug_dir: Path, pts: np.ndarray, centroids: np.ndarray, row_band: np.ndarray, note: str) -> "Path | None":
    """Best-effort scatter colored by attempted row-band, written before a mapping failure is
    raised so there is always something to inspect. Never masks the real error."""
    try:
        plt = _pyplot()
        debug_dir.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(13, 9))
        ax.scatter(pts[:, 1], pts[:, 0], s=6, alpha=0.25, color="#888888", label="field (scene)")
        ax.scatter(centroids[:, 1], centroids[:, 0], s=32, c=row_band, cmap="tab20", label="well centroid")
        for (y, x), rb in zip(centroids, row_band):
            ax.annotate(str(rb), (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7)
        ax.invert_yaxis()
        ax.set_xlabel("stage X (px)")
        ax.set_ylabel("stage Y (px)")
        ax.set_title("Scene -> well mapping FAILED - " + note, fontsize=9)
        ax.legend(loc="upper right", fontsize=8)
        fig.tight_layout()
        path = debug_dir / "scene_to_well_DEBUG.png"
        fig.savefig(path, dpi=150)
        plt.close(fig)
        return path
    except Exception:  # noqa: BLE001
        return None


def write_qc_outputs(out_dir: Path, mapping: SceneWellMapping, field_names: "dict[int, str] | None" = None) -> tuple[Path, Path]:
    """Writes scene_to_well.csv and scene_to_well_qc.png into out_dir; returns both paths."""
    plt = _pyplot()
    out_dir.mkdir(parents=True, exist_ok=True)
    field_names = field_names or {}

    csv_path = out_dir / "scene_to_well.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["Scene", "Field", "Well", "StageY", "StageX", "Source"])
        for (y, x), s in zip(mapping.field_points, mapping.scenes):
            w.writerow([s, field_names.get(s, ""), mapping.scene_to_well[s], int(y), int(x), mapping.source])

    pts, cen = mapping.field_points, mapping.well_centroids
    fig, ax = plt.subplots(figsize=(13, 9))
    ax.scatter(pts[:, 1], pts[:, 0], s=8, alpha=0.35, color="#3b82a0", label="field (scene)")
    ax.scatter(cen[:, 1], cen[:, 0], s=28, color="#6d5bd0", label="well centroid")
    for (y, x), label in zip(cen, mapping.well_labels):
        ax.annotate(label, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("stage X (px)")
    ax.set_ylabel("stage Y (px)")
    ax.set_title(f"Scene -> well mapping QC (source: {mapping.source})")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    png_path = out_dir / "scene_to_well_qc.png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    return csv_path, png_path
