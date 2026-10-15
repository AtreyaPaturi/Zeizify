"""
Command-line interface.

    zeizify                      launch the GUI (same as `zeizify gui`)
    zeizify map     CZI -o OUT   scene -> well mapping + QC only
    zeizify extract CZI -p MAP.csv -o OUT [--yes] [--control AB=WELL:FIELD] [--range AB:CHANNEL=MIN:MAX]
    zeizify info    CZI          print channels, scenes, pixel size, bit depth
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .layout import parse_panel
from .plate import PLATE_FORMATS, normalize_well


def _plate(n: int):
    if n not in PLATE_FORMATS:
        raise argparse.ArgumentTypeError(f"plate must be one of {sorted(PLATE_FORMATS)}")
    return PLATE_FORMATS[n]


def _parse_control(items: list[str]) -> dict[str, tuple[str, str]]:
    out = {}
    for it in items:
        try:
            ab, rest = it.split("=", 1)
            well, fld = rest.split(":", 1)
        except ValueError:
            raise SystemExit(f"--control expects ANTIBODY=WELL:FIELD, got {it!r}")
        out[ab] = (normalize_well(well), fld)
    return out


def _parse_ranges(items: list[str]) -> dict[tuple[str, str], tuple[int, int]]:
    out = {}
    for it in items:
        try:
            key, rng = it.split("=", 1)
            ab, ch = key.rsplit(":", 1)
            lo, hi = (int(v) for v in rng.split(":", 1))
        except ValueError:
            raise SystemExit(f"--range expects ANTIBODY:CHANNEL=MIN:MAX, got {it!r}")
        out[(ab, ch)] = (lo, hi)
    return out


def _parse_colors(items: list[str]) -> dict[str, str]:
    from .layout import COLORS

    out = {}
    for it in items:
        ch, _, color = it.partition("=")
        name = color.strip().capitalize()
        if not ch or name not in COLORS:
            raise SystemExit(f"--color expects CHANNEL=COLOR with COLOR in {', '.join(c.lower() for c in COLORS)}, got {it!r}")
        out[ch.strip()] = name
    return out


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="zeizify", description="Per-well extraction of Zeiss CZI multi-well plates.")
    p.add_argument("--version", action="version", version=f"zeizify {__version__}")
    sub = p.add_subparsers(dest="cmd")
    sub.add_parser("gui", help="launch the desktop app")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("czi", type=Path)
    common.add_argument("--plate", type=lambda s: _plate(int(s)), default=PLATE_FORMATS[96], help="wells per plate (default 96)")
    common.add_argument("--flip-rows", action="store_true", help="geometry fallback only: A1 is at the bottom")
    common.add_argument("--flip-cols", action="store_true", help="geometry fallback only: A1 is at the right")
    common.add_argument("--first-well", default="A1", metavar="WELL",
                        help="geometry fallback only: the imaged well closest to A1, e.g. A2 if imaging started in "
                             "column 2 (default A1)")

    m = sub.add_parser("map", parents=[common], help="scene -> well mapping + QC outputs")
    m.add_argument("-o", "--out", type=Path, required=True)

    e = sub.add_parser("extract", parents=[common], help="full extraction (non-interactive)")
    e.add_argument("-p", "--plate-map", type=Path, required=True)
    e.add_argument("-o", "--out", type=Path, required=True)
    e.add_argument("--yes", action="store_true", help="accept the mapping without the interactive prompt")
    e.add_argument("--control", action="append", default=[], metavar="AB=WELL:FIELD",
                   help="reference field per antibody (default: first control-like treatment)")
    e.add_argument("--range", action="append", default=[], metavar="AB:CH=MIN:MAX",
                   help="hand-entered display range for one antibody/channel")
    e.add_argument("--panel", action="append", default=[], metavar="CH[+CH...]",
                   help="layout panel, e.g. AF488+AF594 for a co-localization merge (repeatable; replaces the defaults)")
    e.add_argument("--color", action="append", default=[], metavar="CH=COLOR",
                   help="pseudocolor for a channel: red, green, blue, magenta, cyan, yellow or gray")
    e.add_argument("--overwrite-png", action="store_true")

    sub.add_parser("info", help="print CZI metadata summary").add_argument("czi", type=Path)
    return p


def main(argv: "list[str] | None" = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd in (None, "gui"):
        from .gui.app import main as gui_main

        gui_main()
        return 0

    from . import pipeline

    if args.cmd == "info":
        from czifile import CziFile

        from .czi_metadata import read_scenes

        with CziFile(str(args.czi)) as czi:
            xml = czi.metadata()
        info = pipeline.read_czi_info(xml)
        scenes = read_scenes(xml)
        print(f"channels:    {info.channel_names}")
        print(f"pixel size:  {info.px_size_um} µm")
        print(f"bit ceiling: {info.bit_ceiling}")
        print(f"scenes:      {len(scenes)} ({sum(bool(s.array_name) for s in scenes)} with ArrayName)")
        return 0

    args.out.mkdir(parents=True, exist_ok=True)
    mapping, png = pipeline.map_scenes(args.czi, args.out, args.plate, args.flip_rows, args.flip_cols,
                                       first_well=args.first_well)
    print(f"QC plot: {png}")
    if args.cmd == "map":
        return 0

    if not args.yes:
        ans = input("Inspect the QC plot. Is the mapping correct? [y/N] ").strip().lower()
        if ans not in ("y", "yes"):
            print("Stopped before processing any wells.")
            return 1

    result = pipeline.extract(args.czi, args.plate_map, args.out, mapping)
    if not result.stats:
        return 1
    from .contrast import default_control_choice

    control = default_control_choice(result.stats)
    control.update(_parse_control(args.control))
    custom = _parse_ranges(args.range)
    panels = [parse_panel(p) for p in args.panel] or None
    colors = _parse_colors(args.color)
    pipeline.finalize(result, args.out, control, custom, args.overwrite_png, panels=panels, colors=colors)
    pipeline.write_manifest(args.out, args.czi, args.plate_map, args.plate, mapping, args.flip_rows, args.flip_cols,
                            result, control, custom, panels, colors)
    return 0


if __name__ == "__main__":
    sys.exit(main())
