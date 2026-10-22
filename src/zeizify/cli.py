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

    from .analysis.params import AnalysisParams

    d = AnalysisParams()
    a = sub.add_parser("analyze", help="per-nucleus spot and colocalization analysis of an extraction")
    a.add_argument("out", type=Path, help="Zeizify output folder (with intensity_report.csv)")
    a.add_argument("--analysis-dir", type=Path, help="where to write results (default OUT/analysis)")
    a.add_argument("--nuclear", default="", help="nuclear channel (default: auto, DAPI/Hoechst)")
    a.add_argument("--channels", default="", help="comma-separated spot channels (default: all non-nuclear)")
    a.add_argument("--method", choices=["auto", "cellpose", "classical"], default=d.nucleus_method)
    a.add_argument("--nucleus-diameter", type=float, default=d.nucleus_diameter_um, metavar="UM")
    a.add_argument("--keep-border-nuclei", action="store_true")
    a.add_argument("--cpu", action="store_true", help="run Cellpose on the CPU")
    a.add_argument("--spot-diameter", type=float, nargs=2, default=[d.spot_min_diameter_um, d.spot_max_diameter_um],
                   metavar=("MIN_UM", "MAX_UM"))
    a.add_argument("--sensitivity", type=float, default=d.spot_sensitivity,
                   help="spot noise floor in MADs (lower = more spots)")
    a.add_argument("--no-declump", action="store_true")
    a.add_argument("--fixed-threshold", action="append", default=[], metavar="CHANNEL[=VALUE]",
                   help="one spot brightness floor for every well in CHANNEL (repeatable); VALUE in top-hat "
                        "intensity, omitted = automatic. For stains whose bodies vanish in some conditions")
    a.add_argument("--dapi-qc", type=float, nargs=2, default=[d.qc_dapi_low, d.qc_dapi_high], metavar=("LOW", "HIGH"),
                   help="flag nuclei whose mean DAPI is outside LOW..HIGH x the field median")
    a.add_argument("--randomizations", type=int, default=d.coloc_randomizations,
                   help="Costes significance randomizations per nucleus (0 = off)")
    a.add_argument("--replicate-column", default="", help="plate-map column defining replicates (default: Replicate, else well)")
    a.add_argument("--treatment-order", default="", metavar="A,B,...",
                   help="comma-separated treatment order in plots and tests (default: control, then plate-map order)")
    a.add_argument("--cell-type-order", default="", metavar="A,B,...",
                   help="comma-separated cell-type order (default: wild type, then plate-map order)")
    a.add_argument("--workers", type=int, default=d.workers)
    a.add_argument("--no-qc", action="store_true", help="skip QC overlay images")

    g = sub.add_parser("figure", help="publication single-cell figure from an analysed output folder")
    g.add_argument("out", type=Path, help="Zeizify output folder, analysed with 'zeizify analyze'")
    g.add_argument("--name", default="single_cells", help="figure name (folder OUT/figures/NAME)")
    g.add_argument("--cell-types", default="", metavar="A,B", help="cell types to show, in order (default: all)")
    g.add_argument("--treatments", default="", metavar="A,B", help="treatments to show, in order (default: all)")
    g.add_argument("--cell-type-labels", default="", metavar="'WT=Wild type; KO1=PML KO'")
    g.add_argument("--treatment-labels", default="", metavar="'ML=ML792; MG=MG132'")
    g.add_argument("--panels", default="", metavar="'DAPI|AF488|AF594|AF488+AF594'",
                   help="panels in order; '+' merges channels (default: nuclear, each signal, all signals merged)")
    g.add_argument("--color", action="append", default=[], metavar="CHANNEL=COLOR")
    g.add_argument("--crop", type=float, nargs=2, default=[30.0, 30.0], metavar=("W_UM", "H_UM"))
    g.add_argument("--rows", choices=["treatment", "panel"], default="treatment",
                   help="treatments as rows (stains across) or stains as rows")
    g.add_argument("--blocks", choices=["side", "stacked"], default="side", help="cell types side by side or stacked")
    g.add_argument("--width", default="journal2",
                   help="journal1 (89 mm), journal1.5 (120 mm), journal2 (183 mm), slide (13.33 in), or a width in mm")
    g.add_argument("--font", default="Arial")
    g.add_argument("--font-size", type=float, default=8.0)
    g.add_argument("--dpi", type=int, default=300)
    g.add_argument("--format", default="png,pdf", help="comma-separated: png, tif, pdf, svg")
    g.add_argument("--contrast", choices=["shared", "panel"], default="shared",
                   help="one display range per channel for all conditions, or per panel")
    g.add_argument("--saturated", type=float, default=0.3, help="%% of nuclear pixels allowed to saturate")
    g.add_argument("--max", action="append", default=[], metavar="CHANNEL=VALUE", help="fixed display maximum")
    g.add_argument("--outline", choices=["none", "nuclear", "all"], default="nuclear")
    g.add_argument("--scale-bar", type=float, default=5.0, metavar="UM", help="0 = none")
    g.add_argument("--pick", action="append", default=[], metavar="CELLTYPE,TREATMENT=WELL,FIELD,NUCLEUS",
                   help="use this nucleus instead of the automatic pick (repeatable)")
    return p


def _figure(args) -> int:
    from . import cellfigure as cf

    cd = cf.load_cells(args.out)
    cts = _csv(args.cell_types) or list(dict.fromkeys(c for c, _ in cd.conditions()))
    trs = _csv(args.treatments) or list(dict.fromkeys(t for _, t in cd.conditions()))
    widths = {"journal1": 89 / 25.4, "journal1.5": 120 / 25.4, "journal2": 183 / 25.4, "slide": 13.333}
    width = widths.get(args.width.lower()) or float(args.width) / 25.4
    colors = _parse_colors(args.color)
    if cd.nuclear and cd.nuclear not in colors:
        colors[cd.nuclear] = "Gray"
    opts = cf.FigureOptions(
        panels=cf.parse_order(args.panels) if args.panels else cf.default_figure_panels(cd.ds.channels),
        colors=colors, cell_types=cts, treatments=trs, cell_type_labels=cf.parse_labels(args.cell_type_labels),
        treatment_labels=cf.parse_labels(args.treatment_labels), crop_w_um=args.crop[0], crop_h_um=args.crop[1],
        rows=args.rows, blocks=args.blocks, width_in=width, dpi=args.dpi, font=args.font, font_size=args.font_size,
        contrast=args.contrast, saturated_pct=args.saturated, max_override={k: v for k, v in _fixed(args.max).items() if v},
        outline=args.outline, scale_bar_um=args.scale_bar, scale_bar="last" if args.scale_bar > 0 else "none",
        formats=_csv(args.format))
    conds = [(c, t) for c in cts for t in trs if (c, t) in set(cd.conditions())]
    manual = {}
    for it in args.pick:
        cond, _, where = it.partition("=")
        ct, _, t = cond.partition(",")
        well, fld, nid = [x.strip() for x in where.split(",")]
        manual[(ct.strip(), t.strip())] = cf.CellChoice(ct.strip(), t.strip(), well, fld, int(nid), "manual")
    ranked = cf.auto_choices(cd, [k for k in conds if k not in manual], opts.crop_w_um, opts.crop_h_um, n=5,
                             progress=lambda i, n, item: print(f"  ranked {i}/{n}: {item}"))
    choices = {}
    for k in conds:
        if k in manual:
            choices[k] = manual[k]
        elif len(ranked.get(k, [])):
            choices[k] = cf.choice_from_row(ranked[k].iloc[0])
            print(f"{k[0]} · {k[1]}: {cf.describe(ranked[k].iloc[0], cd)}")
        else:
            print(f"{k[0]} · {k[1]}: no suitable cell - shown as n/a (use --pick)")
    res = cf.export(cd, choices, opts, args.name)
    print("Wrote", ", ".join(str(p) for p in res["files"]))
    print(f"Cells, legend and settings: {res['folder']}")
    return 0


def _csv(text: str) -> list[str]:
    return [t.strip() for t in text.split(",") if t.strip()]


def _fixed(items: list[str]) -> dict[str, float]:
    out = {}
    for it in items:
        ch, _, val = it.partition("=")
        out[ch.strip()] = float(val) if val.strip() else 0.0
    return out


def _analyze(args) -> int:
    from .analysis.params import AnalysisParams
    from .analysis.run import run_analysis

    p = AnalysisParams(
        nuclear_channel=args.nuclear, spot_channels=_csv(args.channels),
        nucleus_method=args.method, nucleus_diameter_um=args.nucleus_diameter,
        exclude_border_nuclei=not args.keep_border_nuclei, use_gpu=not args.cpu,
        spot_min_diameter_um=args.spot_diameter[0], spot_max_diameter_um=args.spot_diameter[1],
        spot_sensitivity=args.sensitivity, spot_declump=not args.no_declump,
        qc_dapi_low=args.dapi_qc[0], qc_dapi_high=args.dapi_qc[1], coloc_randomizations=args.randomizations,
        replicate_column=args.replicate_column, workers=args.workers, qc_overlays=not args.no_qc,
        treatment_order=_csv(args.treatment_order), cell_type_order=_csv(args.cell_type_order),
        spot_fixed_threshold=_fixed(args.fixed_threshold))
    last = {"stage": None}

    def progress(stage, done, total, item):
        if stage != last["stage"]:
            last["stage"] = stage
            print(f"== {stage} ==")

    res = run_analysis(args.out, p, args.analysis_dir, progress=progress)
    for k, e in res.errors:
        print(f"  error in {k}: {e}")
    print(f"Results: {res.analysis_dir}")
    return 1 if res.errors and not len(res.nuclei) else 0


def main(argv: "list[str] | None" = None) -> int:
    args = build_parser().parse_args(argv)
    if args.cmd in (None, "gui"):
        from .gui.app import main as gui_main

        gui_main()
        return 0

    if args.cmd == "analyze":
        return _analyze(args)
    if args.cmd == "figure":
        return _figure(args)

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
