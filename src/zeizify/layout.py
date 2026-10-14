"""
Figure-style layout PNGs: one image per field with solo and DAPI-merged panels side by
side, contrast set by the antibody's reference range. Built by reading back the small
per-channel TIFFs - the CZI is never re-decoded and the TIFFs are never rewritten.
"""

from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont

from .contrast import ControlChoice, CustomRanges, build_control_range
from .extraction import FieldStat
from .plate import slug

RGB = tuple[float, float, float]

# Pseudocolor by fluorophore name (substring match, case-insensitive). First hit wins.
CHANNEL_COLOR_RULES: list[tuple[tuple[str, ...], RGB]] = [
    (("dapi", "hoechst", "405", "bfp"), (0.0, 0.0, 1.0)),
    (("488", "gfp", "fitc", "alexa488"), (0.0, 1.0, 0.0)),
    (("647", "cy5", "far red", "farred", "640", "670"), (1.0, 0.0, 1.0)),
    (("594", "568", "555", "561", "cy3", "tritc", "rfp", "mcherry", "tdtomato", "texas"), (1.0, 0.0, 0.0)),
]
NUCLEAR_KEYWORDS = ("dapi", "hoechst")

LAYOUT_GAP_PX = 24
LAYOUT_GAP_COLOR = (255, 255, 255)
LAYOUT_LABEL_FONT_SIZE = 48


def channel_color(name: str) -> RGB:
    low = name.lower()
    for keys, rgb in CHANNEL_COLOR_RULES:
        if any(k in low for k in keys):
            return rgb
    return (1.0, 1.0, 1.0)  # unknown -> grayscale


def default_panels(channel_names: list[str]) -> list[tuple[str, ...]]:
    """Every non-nuclear channel solo, then the nuclear stain solo, then each non-nuclear
    channel merged with it. For AF488/AF594/DAPI this gives
    AF488 | AF594 | DAPI | AF488+DAPI | AF594+DAPI."""
    nuclear = next((c for c in channel_names if any(k in c.lower() for k in NUCLEAR_KEYWORDS)), None)
    others = [c for c in channel_names if c != nuclear]
    panels: list[tuple[str, ...]] = [(c,) for c in others]
    if nuclear:
        panels.append((nuclear,))
        panels += [(c, nuclear) for c in others]
    return panels


def scale_to_uint8(arr: np.ndarray, lo: float, hi: float) -> np.ndarray:
    if hi <= lo:
        hi = lo + 1
    x = (arr.astype(np.float32) - lo) / (hi - lo)
    return (np.clip(x, 0.0, 1.0) * 255).astype(np.uint8)


def compose_rgb(scaled: dict[str, np.ndarray], spec: tuple[str, ...]) -> np.ndarray:
    """Additive pseudocolor blend of the channels in spec, clipped to 8 bit."""
    acc = np.zeros((*scaled[spec[0]].shape, 3), dtype=np.float32)
    for cname in spec:
        acc += scaled[cname][..., None].astype(np.float32) * np.array(channel_color(cname), dtype=np.float32)
    return np.clip(acc, 0, 255).astype(np.uint8)


def _label_panel(img: Image.Image, text: str) -> None:
    draw = ImageDraw.Draw(img)
    size = max(12, min(LAYOUT_LABEL_FONT_SIZE, img.width // 12))
    try:
        font = ImageFont.load_default(size=size)
    except TypeError:  # Pillow < 10.1
        font = ImageFont.load_default()
    pad = max(4, size // 5)
    x0, y0, x1, y1 = draw.textbbox((0, 0), text, font=font)
    draw.rectangle([0, 0, (x1 - x0) + 2 * pad, (y1 - y0) + 2 * pad], fill=(0, 0, 0))
    draw.text((pad, pad - y0), text, fill=(255, 255, 255), font=font)


def generate_png_layouts(
    stats: list[FieldStat],
    control: ControlChoice,
    out_dir: Path,
    custom: "CustomRanges | None" = None,
    panels: "list[tuple[str, ...]] | None" = None,
    overwrite: bool = False,
) -> int:
    """One PNG per field at OUT_DIR/<Antibody>/png/<Label>_<Well>_<Field>_layout.png.
    Existing PNGs are kept unless overwrite=True (e.g. after changing the contrast)."""
    rng = build_control_range(stats, control, custom)
    by_field: dict[tuple, dict[str, FieldStat]] = defaultdict(dict)
    for s in stats:
        by_field[(s.well, s.field, s.antibody, s.label)][s.channel_name] = s

    written = 0
    for (well, field, antibody, label), channels in by_field.items():
        folder = out_dir / (slug(antibody) or "unassigned") / "png"
        folder.mkdir(parents=True, exist_ok=True)
        out_path = folder / f"{label}_{well}_{field}_layout.png"
        if out_path.exists() and not overwrite:
            continue

        scaled = {}
        for cname, s in channels.items():
            lo, hi = rng.get((antibody, cname), (s.min, s.max))
            scaled[cname] = scale_to_uint8(tifffile.imread(s.path), lo, hi)

        images = []
        for spec in panels or default_panels(list(channels)):
            if all(c in scaled for c in spec):
                img = Image.fromarray(compose_rgb(scaled, spec))
                _label_panel(img, "+".join(spec))
                images.append(img)
        if not images:
            continue

        pw, ph = images[0].size
        canvas = Image.new("RGB", (pw * len(images) + LAYOUT_GAP_PX * (len(images) - 1), ph), LAYOUT_GAP_COLOR)
        for k, img in enumerate(images):
            canvas.paste(img, (k * (pw + LAYOUT_GAP_PX), 0))
        canvas.save(out_path)
        written += 1
    return written
