"""
Figure-style layout PNGs: one image per field with the chosen panels side by side, contrast
set by the antibody panel's reference range. Built by reading back the small per-channel
TIFFs - the CZI is never re-decoded and the TIFFs are never rewritten.

Panels are tuples of channel names: ("AF488",) is a single channel, ("AF488", "AF594") a
co-localization merge, ("AF488", "AF594", "DAPI") a merge with the nuclear stain. Users
choose the panels and each channel's pseudocolor; defaults are inferred from channel names.
"""

from __future__ import annotations

import threading
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Callable

import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont

from .contrast import ControlChoice, CustomRanges, build_control_range
from .extraction import FieldStat
from .plate import slug

RGB = tuple[float, float, float]
Panel = tuple[str, ...]

COLORS: dict[str, RGB] = {
    "Red": (1.0, 0.0, 0.0),
    "Green": (0.0, 1.0, 0.0),
    "Blue": (0.0, 0.0, 1.0),
    "Magenta": (1.0, 0.0, 1.0),
    "Cyan": (0.0, 1.0, 1.0),
    "Yellow": (1.0, 1.0, 0.0),
    "Gray": (1.0, 1.0, 1.0),
}

# Default pseudocolor by fluorophore name (substring match, case-insensitive). First hit wins.
CHANNEL_COLOR_RULES: list[tuple[tuple[str, ...], str]] = [
    (("dapi", "hoechst", "405", "bfp"), "Blue"),
    (("488", "gfp", "fitc", "alexa488"), "Green"),
    (("647", "cy5", "far red", "farred", "640", "670"), "Magenta"),
    (("594", "568", "555", "561", "546", "cy3", "tritc", "rfp", "mcherry", "tdtomato", "texas"), "Red"),
]
NUCLEAR_KEYWORDS = ("dapi", "hoechst")

LAYOUT_GAP_PX = 24
LAYOUT_GAP_COLOR = (255, 255, 255)
LAYOUT_LABEL_FONT_SIZE = 48


def default_color_name(name: str) -> str:
    low = name.lower()
    for keys, color in CHANNEL_COLOR_RULES:
        if any(k in low for k in keys):
            return color
    return "Gray"


def channel_color(name: str, colors: "dict[str, str] | None" = None) -> RGB:
    """RGB weights for a channel: the user's choice if given, else the name-based default."""
    return COLORS[(colors or {}).get(name) or default_color_name(name)]


def nuclear_channel(channel_names: list[str]) -> "str | None":
    return next((c for c in channel_names if any(k in c.lower() for k in NUCLEAR_KEYWORDS)), None)


def candidate_panels(channel_names: list[str]) -> list[Panel]:
    """Every panel worth offering, in display order: each channel alone, each signal channel
    merged with the nuclear stain, every co-localization combination of signal channels, and
    each such combination merged with the nuclear stain."""
    nuc = nuclear_channel(channel_names)
    sig = [c for c in channel_names if c != nuc]
    out: list[Panel] = [(c,) for c in sig] + ([(nuc,)] if nuc else [])
    if nuc:
        out += [(c, nuc) for c in sig]
    for k in range(2, len(sig) + 1):
        for combo in combinations(sig, k):
            out.append(combo)
            if nuc:
                out.append((*combo, nuc))
    return out


def default_panels(channel_names: list[str]) -> list[Panel]:
    """Each signal channel alone, the nuclear stain, each signal + nuclear, and - when there are
    two or more signal channels - all signal channels merged, with and without the nuclear stain.
    For AF488/AF594/DAPI: AF488 | AF594 | DAPI | AF488+DAPI | AF594+DAPI | AF488+AF594 | AF488+AF594+DAPI."""
    nuc = nuclear_channel(channel_names)
    sig = [c for c in channel_names if c != nuc]
    panels: list[Panel] = [(c,) for c in sig]
    if nuc:
        panels.append((nuc,))
        panels += [(c, nuc) for c in sig]
    if len(sig) >= 2:
        panels.append(tuple(sig))
        if nuc:
            panels.append((*sig, nuc))
    return panels


def parse_panel(text: str) -> Panel:
    """'AF488+AF594' -> ('AF488', 'AF594')."""
    return tuple(p.strip() for p in text.split("+") if p.strip())


def scale_to_uint8(arr: np.ndarray, lo: float, hi: float) -> np.ndarray:
    if hi <= lo:
        hi = lo + 1
    x = (arr.astype(np.float32) - lo) / (hi - lo)
    return (np.clip(x, 0.0, 1.0) * 255).astype(np.uint8)


def compose_rgb(scaled: dict[str, np.ndarray], spec: Panel, colors: "dict[str, str] | None" = None) -> np.ndarray:
    """Additive pseudocolor blend of the channels in spec, clipped to 8 bit (so overlapping
    green + magenta signal shows as white)."""
    acc = np.zeros((*scaled[spec[0]].shape, 3), dtype=np.float32)
    for cname in spec:
        acc += scaled[cname][..., None].astype(np.float32) * np.array(channel_color(cname, colors), dtype=np.float32)
    return np.clip(acc, 0, 255).astype(np.uint8)


def panel_caption(spec: Panel, stains: dict[str, str]) -> str:
    """'CD7 (AF488)' for a single stained channel; 'CD7 + HLA + DAPI' for merges."""
    if len(spec) == 1:
        c = spec[0]
        return f"{stains[c]} ({c})" if stains.get(c) else c
    return " + ".join(stains.get(c) or c for c in spec)


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


def stains_by_channel(stats: list[FieldStat]) -> dict[str, dict[str, str]]:
    """{antibody panel: {channel: antibody}} - the most common stain per channel within a panel."""
    counts: dict[str, dict[str, Counter]] = defaultdict(lambda: defaultdict(Counter))
    for s in stats:
        if s.stain:
            counts[s.antibody][s.channel_name][s.stain] += 1
    return {ab: {ch: c.most_common(1)[0][0] for ch, c in chs.items()} for ab, chs in counts.items()}


def generate_png_layouts(
    stats: list[FieldStat],
    control: ControlChoice,
    out_dir: Path,
    custom: "CustomRanges | None" = None,
    panels: "list[Panel] | None" = None,
    overwrite: bool = False,
    colors: "dict[str, str] | None" = None,
    progress: "Callable[[int, int, str], None] | None" = None,
    log: "Callable[[str], None] | None" = None,
    cancel: "threading.Event | None" = None,
) -> int:
    """One PNG per field at OUT_DIR/<Antibody>/png/<Label>_<Well>_<Field>_layout.png.
    Panels a field lacks a channel for are skipped. Existing PNGs are kept unless overwrite=True
    (e.g. after changing the contrast, panels or colors).

    progress(done, total, name) is called after every field (written or kept), log(text) once per
    field with what happened, and a set cancel event stops after the current field. Returns the
    number of PNGs written."""
    rng = build_control_range(stats, control, custom)
    by_field: dict[tuple, dict[str, FieldStat]] = defaultdict(dict)
    for s in stats:
        by_field[(s.well, s.field, s.antibody, s.label)][s.channel_name] = s

    total = len(by_field)
    written = 0
    for i, ((well, field, antibody, label), channels) in enumerate(by_field.items(), start=1):
        if cancel is not None and cancel.is_set():
            if log:
                log(f"Stopped after {i - 1} of {total} fields - finished PNGs are kept.")
            break
        name = f"{label}_{well}_{field}_layout.png"
        result = _write_layout(out_dir, well, field, antibody, label, channels, rng, panels, overwrite, colors)
        written += result == "written"
        if log:
            log({"written": f"wrote   {name}", "kept": f"kept    {name} (already exists)",
                 "empty": f"skipped {name} (none of the chosen panels apply)"}[result])
        if progress:
            progress(i, total, f"{well} {field} · {antibody}")
    return written


def _write_layout(out_dir, well, field, antibody, label, channels, rng, panels, overwrite, colors) -> str:
    """Write one field's layout PNG. Returns "written", "kept" (exists, no overwrite) or "empty"."""
    folder = out_dir / (slug(antibody) or "unassigned") / "png"
    folder.mkdir(parents=True, exist_ok=True)
    out_path = folder / f"{label}_{well}_{field}_layout.png"
    if out_path.exists() and not overwrite:
        return "kept"

    scaled = {}
    for cname, s in channels.items():
        lo, hi = rng.get((antibody, cname), (s.min, s.max))
        scaled[cname] = scale_to_uint8(tifffile.imread(s.path), lo, hi)
    stains = {c: s.stain for c, s in channels.items() if s.stain}

    images = []
    for spec in panels or default_panels(list(channels)):
        if spec and all(c in scaled for c in spec):
            img = Image.fromarray(compose_rgb(scaled, spec, colors))
            _label_panel(img, panel_caption(spec, stains))
            images.append(img)
    if not images:
        return "empty"

    pw, ph = images[0].size
    canvas = Image.new("RGB", (pw * len(images) + LAYOUT_GAP_PX * (len(images) - 1), ph), LAYOUT_GAP_COLOR)
    for k, img in enumerate(images):
        canvas.paste(img, (k * (pw + LAYOUT_GAP_PX), 0))
    canvas.save(out_path)
    return "written"
