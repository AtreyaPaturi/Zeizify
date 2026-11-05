"""
Run the analysis over every field of a Zeizify output folder, in parallel.

    stage 1  nuclei    Cellpose on the GPU in this process (one model, no GPU contention), or the
                       classical method on all CPU cores. Masks are cached in analysis/masks and reused
                       while the nucleus settings are unchanged - re-tuning spot settings is fast.
    stage 2  measure   spots, per-nucleus measurements and colocalization, on all CPU cores.
    stage 3  stats     tables, replicate-level statistics, SuperPlots, manifest.

Outputs (OUT_DIR/analysis/):
    nuclei.csv  spots.csv  colocalization.csv  fields.csv
    replicate_summary.csv  condition_summary.csv  stats_omnibus.csv  stats_pairwise.csv
    superplots/  qc/  masks/  analysis_manifest.json
"""

from __future__ import annotations

import hashlib
import json
import platform
import threading
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import tifffile

from .. import __version__
from .io import Dataset, FieldImages, load_dataset
from .measure import FieldResult, measure_field
from .params import AnalysisParams
from .segment import (SEPARABILITY_MIN, auto_spot_threshold, detect_spots, resolve_method, segment_nuclei,
                      spot_brightness)

Log = Callable[[str], None]
Progress = Callable[[str, int, int, str], None]  # (stage, done, total, item)


@dataclass
class AnalysisOutput:
    analysis_dir: Path
    nuclei: pd.DataFrame
    spots: pd.DataFrame
    coloc: pd.DataFrame
    condition_summary: pd.DataFrame
    omnibus: pd.DataFrame
    pairwise: pd.DataFrame
    overview_plots: list[Path]
    errors: list[tuple[str, str]]
    method: str


def _nuclei_signature(params: AnalysisParams, method: str, nuclear: str) -> str:
    keys = (method, nuclear, params.nucleus_diameter_um, params.nucleus_min_area_frac, params.nucleus_max_area_frac,
            params.exclude_border_nuclei)
    return hashlib.sha1(repr(keys).encode()).hexdigest()[:12]


def _mask_path(masks: Path, fi: FieldImages) -> Path:
    return masks / f"{fi.key}_nuclei.tif"


def _cells_path(masks: Path, fi: FieldImages) -> Path:
    return masks / f"{fi.key}_cells.tif"


def _cells_signature(params: AnalysisParams, method: str, nuclei_sig: str) -> str:
    keys = (nuclei_sig, method, params.cell_marker, params.cell_diameter_um, params.exclude_border_cells)
    return hashlib.sha1(repr(keys).encode()).hexdigest()[:12]


def _cells_task(fi: FieldImages, marker: str, nuclear: str, nuclei_path: Path, params: AnalysisParams, px_um: float,
                method: str, out: Path) -> str:
    from .cells import segment_cells

    nuclei = tifffile.imread(nuclei_path).astype(np.int32)
    lab = segment_cells(tifffile.imread(fi.channels[marker]), tifffile.imread(fi.channels[nuclear]), nuclei, params,
                        px_um, method)
    tifffile.imwrite(out, lab.astype(np.uint16 if lab.max() < 65535 else np.uint32), compression="zlib")
    return fi.key


def _segment_task(fi: FieldImages, nuclear: str, params: AnalysisParams, px_um: float, method: str, out: Path) -> str:
    lab = segment_nuclei(tifffile.imread(fi.channels[nuclear]), params, px_um, method)
    tifffile.imwrite(out, lab.astype(np.uint16 if lab.max() < 65535 else np.uint32), compression="zlib")
    return fi.key


def _brightness_task(fi: FieldImages, mask: Path, channels: list[str], params: AnalysisParams,
                     px_um: float) -> dict[str, np.ndarray]:
    """Mean top-hat brightness of the per-nucleus spot candidates in one field (for the automatic
    experiment-wide floor)."""
    import warnings

    warnings.filterwarnings("ignore")
    nuclei = tifffile.imread(mask).astype(np.int32)
    out = {}
    for c in channels:
        if c in fi.channels:
            out[c] = spot_brightness(*detect_spots(tifffile.imread(fi.channels[c]), nuclei, params, px_um))
    return out


def resolve_spot_floors(params: AnalysisParams, channels: list[str], candidates: dict) -> dict[str, dict]:
    """Per fixed-threshold channel: {"floor", "mode", "separability"}; automatic floors come from the
    pooled `candidates` (channel -> array of spot brightness)."""
    out = {}
    for c, v in params.spot_fixed_threshold.items():
        if c not in channels:
            continue
        if v and v > 0:
            out[c] = {"floor": float(v), "mode": "fixed (manual)", "separability": float("nan")}
        else:
            t, eta = auto_spot_threshold(np.asarray(candidates.get(c, []), dtype=float))
            out[c] = {"floor": t, "mode": "fixed (automatic)", "separability": eta}
    return out


def _actin_hist_task(fi: FieldImages, mask: Path, cells_path: "Path | None", channel: str, px_um: float,
                     ceiling: float) -> np.ndarray:
    import warnings

    warnings.filterwarnings("ignore")
    from .actin import reference_histogram

    return reference_histogram(fi, mask, cells_path, channel, px_um, ceiling)


def resolve_actin_thresholds(ds: Dataset, ready: list, channels: list[str], params: AnalysisParams, masks: Path,
                             markers: dict, log: Log) -> tuple[dict, pd.DataFrame]:
    """Bright F-actin cut-off per (staining panel, channel): the `actin_bright_percentile` of the reference
    condition's cell pixels (background-subtracted), or the fixed `actin_bright_threshold`. Returns
    ({field key: {channel: cut-off}}, table for actin_thresholds.csv)."""
    from .actin import is_actin, threshold_from_histogram
    from .stats import default_orders, resolve_reference

    meta = pd.DataFrame([{"cell_type": f.cell_type or "", "treatment": f.treatment or "", "well": f.well}
                         for f in ds.fields])
    ref_ct, ref_t = resolve_reference(default_orders(meta, params.treatment_order, params.cell_type_order),
                                      params.reference_cell_type, params.reference_treatment)
    ceiling = float(ds.bit_ceiling or 65535)
    per_field: dict[str, dict] = {}
    rows = []
    for c in channels:
        stained = [f for f in ready if is_actin(f.stains.get(c, ""))]
        for ab in dict.fromkeys(f.antibody for f in stained):
            panel = [f for f in stained if f.antibody == ab]
            ref = [f for f in panel if (f.cell_type or "") == ref_ct and (f.treatment or "") == ref_t]
            note = ""
            if params.actin_bright_threshold and params.actin_bright_threshold > 0:
                thr, mode = float(params.actin_bright_threshold), "fixed"
            else:
                if not ref:
                    ref, note = panel, "reference condition not in this panel - all its fields pooled"
                hists = []
                with ProcessPoolExecutor(max_workers=params.workers) as ex:
                    futs = [ex.submit(_actin_hist_task, f, _mask_path(masks, f),
                                      _cells_path(masks, f) if markers.get(f.key) else None, c, ds.px_um, ceiling)
                            for f in ref]
                    for fut in futs:
                        try:
                            hists.append(fut.result())
                        except Exception:  # noqa: BLE001 - a bad field is reported by the measure step
                            pass
                thr = (threshold_from_histogram(np.sum(hists, axis=0), ceiling, params.actin_bright_percentile)
                       if hists else float("nan"))
                mode = f"P{params.actin_bright_percentile:g} of {ref_ct} · {ref_t} cell pixels"
            for f in panel:
                per_field.setdefault(f.key, {})[c] = thr
            rows.append({"antibody": ab, "channel": c, "stain": panel[0].stains.get(c, ""), "threshold": thr,
                         "mode": mode, "reference_fields": len(ref), "note": note})
            log(f"  {c} ({panel[0].stains.get(c, '')}, {ab}): bright F-actin above {thr:.0f} ({mode})"
                + (f" - {note}" if note else ""))
    return per_field, pd.DataFrame(rows)


def _measure_task(fi: FieldImages, mask: Path, nuclear: str, spot_channels: list[str], params: AnalysisParams,
                  px_um: float, bits: "int | None", qc_png: "Path | None",
                  floors: "dict[str, float] | None" = None, cells_path: "Path | None" = None,
                  marker: "str | None" = None, actin: "dict[str, float] | None" = None) -> FieldResult:
    import warnings

    warnings.filterwarnings("ignore")
    try:
        nuclei = tifffile.imread(mask).astype(np.int32)
        cells = tifffile.imread(cells_path).astype(np.int32) if cells_path is not None and cells_path.exists() else None
        return measure_field(fi, nuclei, nuclear, spot_channels, params, px_um, bits, qc_png, floors, cells,
                             marker if cells is not None else None, actin)
    except Exception as exc:  # noqa: BLE001 - one bad field never aborts the run
        return FieldResult(fi.key, error=f"{exc}\n{traceback.format_exc(limit=3)}")


def cellpose_device(params: AnalysisParams) -> str:
    """Human-readable device Cellpose will use, with a hint when a GPU was requested but unavailable."""
    if not params.use_gpu:
        return "CPU (GPU disabled)"
    try:
        import torch

        if torch.cuda.is_available():
            return f"GPU ({torch.cuda.get_device_name(0)})"
        return ("CPU - WARNING: PyTorch has no CUDA support in this environment; install the CUDA build "
                "(see README, 'GPU') for ~10x faster segmentation")
    except Exception:  # noqa: BLE001
        return "CPU (PyTorch unavailable)"


def replicate_unit(params: AnalysisParams, ds: Dataset) -> str:
    rc = params.replicate_column.lower()
    if rc == "field":
        return "field (technical replicate - fields of one well are not independent)"
    if rc == "well" or (not rc and all(f.replicate == f.well for f in ds.fields)):
        return "well"
    return f"plate-map column {params.replicate_column or 'Replicate'}"


def _versions() -> dict:
    out = {"zeizify": __version__, "python": platform.python_version()}
    for pkg in ("numpy", "scipy", "scikit-image", "pandas", "tifffile", "cellpose", "torch"):
        try:
            out[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return out


def run_analysis(out_dir: "str | Path", params: "AnalysisParams | None" = None, analysis_dir: "Path | None" = None,
                 log: Log = print, progress: "Progress | None" = None, cancel: "threading.Event | None" = None,
                 fields: "list[str] | None" = None) -> AnalysisOutput:
    params = params or AnalysisParams()
    t0 = time.monotonic()
    ds: Dataset = load_dataset(out_dir, params.replicate_column)
    if not ds.fields:
        raise ValueError("No fields found in the intensity report.")
    if not ds.px_um:
        raise ValueError("Pixel size unknown (not in run_manifest.json or the TIFFs) - cannot convert µm settings.")
    nuclear = params.nuclear_channel or ds.nuclear
    if not nuclear:
        raise ValueError(f"No nuclear channel found among {ds.channels}; set it explicitly.")
    spot_channels = params.spot_channels or ds.signal_channels
    todo = [f for f in ds.fields if not fields or f.key in fields]
    method = resolve_method(params)
    adir = Path(analysis_dir or Path(out_dir) / "analysis")
    masks = adir / "masks"
    masks.mkdir(parents=True, exist_ok=True)
    stop = cancel or threading.Event()
    prog = progress or (lambda *a: None)

    log(f"{len(todo)} field(s) · channels: nuclear {nuclear}, spots {', '.join(spot_channels)} · "
        f"{ds.px_um:.4f} µm/px · nuclei: {method} · {params.workers} worker(s)")
    if method == "cellpose":
        log(f"Cellpose device: {cellpose_device(params)}")
    timings: dict[str, float] = {}
    t_stage = time.monotonic()

    # ---- stage 1: nuclei (cached) ----
    sig = _nuclei_signature(params, method, nuclear)
    sig_file = masks / "signature.json"
    old = json.loads(sig_file.read_text()) if sig_file.exists() else {}
    if old.get("signature") != sig:
        for p in masks.glob("*_nuclei.tif"):
            p.unlink()
        sig_file.write_text(json.dumps({"signature": sig, "method": method, "params": params.to_dict()}, indent=2))
    need = [f for f in todo if not _mask_path(masks, f).exists()]
    log(f"Nuclei: {len(todo) - len(need)} cached, {len(need)} to segment.")
    done = len(todo) - len(need)
    prog("nuclei", done, len(todo), "")
    errors: list[tuple[str, str]] = []
    if need and method == "cellpose":
        for f in need:
            if stop.is_set():
                break
            try:
                _segment_task(f, nuclear, params, ds.px_um, method, _mask_path(masks, f))
            except Exception as exc:  # noqa: BLE001
                errors.append((f.key, f"nuclei: {exc}"))
            done += 1
            prog("nuclei", done, len(todo), f"{f.well} {f.field}")
    elif need:
        with ProcessPoolExecutor(max_workers=params.workers) as ex:
            futs = {ex.submit(_segment_task, f, nuclear, params, ds.px_um, method, _mask_path(masks, f)): f for f in need}
            for fut in as_completed(futs):
                f = futs[fut]
                if stop.is_set():
                    for x in futs:
                        x.cancel()
                try:
                    fut.result()
                except Exception as exc:  # noqa: BLE001
                    errors.append((f.key, f"nuclei: {exc}"))
                done += 1
                prog("nuclei", done, len(todo), f"{f.well} {f.field}")
    if stop.is_set():
        raise InterruptedError("Stopped during nuclear segmentation - finished masks are cached.")
    timings["nuclei_s"] = round(time.monotonic() - t_stage, 1)
    log(f"Nuclei done in {timings['nuclei_s']} s.")
    t_stage = time.monotonic()

    # ---- stage 1b: whole cells (fields with a cell-boundary marker, e.g. F-actin), cached ----
    from .cells import cell_marker_for, resolve_cell_method

    markers = {f.key: cell_marker_for(f.stains, f.channels, params) for f in todo}
    cell_fields = [f for f in todo if markers[f.key] and _mask_path(masks, f).exists()]
    cell_method = resolve_cell_method(params) if cell_fields else ""
    if "cytoplasm" in (params.spot_compartments or []) and len(cell_fields) < len(todo):
        log(f"Note: cytoplasmic spots need cell boundaries - {len(todo) - len(cell_fields)} of {len(todo)} field(s) "
            "have no cell-boundary stain (e.g. F-actin), so only their nuclei are quantified.")
    if cell_fields:
        csig = _cells_signature(params, cell_method, sig)
        csig_file = masks / "cells_signature.json"
        cold = json.loads(csig_file.read_text()) if csig_file.exists() else {}
        if cold.get("signature") != csig:
            for p in masks.glob("*_cells.tif"):
                p.unlink()
            csig_file.write_text(json.dumps({"signature": csig, "method": cell_method}, indent=2))
        need_c = [f for f in cell_fields if not _cells_path(masks, f).exists()]
        used = sorted({f"{markers[f.key]} ({f.stains.get(markers[f.key], '')})" for f in cell_fields})
        log(f"Cells: {len(cell_fields)} field(s) with a cell-boundary marker - {', '.join(used)} - method "
            f"{cell_method}; {len(cell_fields) - len(need_c)} cached, {len(need_c)} to segment.")
        args = lambda f: (f, markers[f.key], nuclear, _mask_path(masks, f), params, ds.px_um, cell_method,  # noqa: E731
                          _cells_path(masks, f))
        if need_c and cell_method == "cellpose":
            for i, f in enumerate(need_c, start=1):
                if stop.is_set():
                    break
                try:
                    _cells_task(*args(f))
                except Exception as exc:  # noqa: BLE001
                    errors.append((f.key, f"cells: {exc}"))
                prog("nuclei", i, len(need_c), f"cells {f.well} {f.field}")
        elif need_c:
            with ProcessPoolExecutor(max_workers=params.workers) as ex:
                futs = {ex.submit(_cells_task, *args(f)): f for f in need_c}
                for i, fut in enumerate(as_completed(futs), start=1):
                    f = futs[fut]
                    if stop.is_set():
                        for x in futs:
                            x.cancel()
                    try:
                        fut.result()
                    except Exception as exc:  # noqa: BLE001
                        errors.append((f.key, f"cells: {exc}"))
                    prog("nuclei", i, len(need_c), f"cells {f.well} {f.field}")
        if stop.is_set():
            raise InterruptedError("Stopped during cell segmentation - finished masks are cached.")
        timings["cells_s"] = round(time.monotonic() - t_stage, 1)
        log(f"Cells done in {timings['cells_s']} s.")
        t_stage = time.monotonic()

    # ---- stage 2: measure (parallel) ----
    qc_dir = adir / "qc"
    results: list[FieldResult] = []
    ready = [f for f in todo if _mask_path(masks, f).exists()]

    # experiment-wide spot floors: pool the per-nucleus spot candidates of every field first
    speckle_ch = [c for c in spot_channels if "speckles" in params.quantify_for(c)]
    signal_ch = [c for c in spot_channels if "signal" in params.quantify_for(c)]
    do_spots, do_signal = bool(speckle_ch), bool(signal_ch)
    log("Quantifying: " + "; ".join(
        f"{c} " + " + ".join(x for x, on in (("speckles", c in speckle_ch), ("overall signal", c in signal_ch)) if on)
        for c in spot_channels) + ".")
    auto = [c for c, v in params.spot_fixed_threshold.items() if c in speckle_ch and not (v and v > 0)]
    candidates: dict[str, list] = {c: [] for c in auto}
    tagged: list[pd.DataFrame] = []  # the same, with treatment / cell type, for the check plot
    if auto:
        log(f"Fixed spot threshold (automatic) for {', '.join(auto)}: pooling spot brightness over all fields...")
        with ProcessPoolExecutor(max_workers=params.workers) as ex:
            futs = {ex.submit(_brightness_task, f, _mask_path(masks, f), auto, params, ds.px_um): f for f in ready}
            for i, fut in enumerate(as_completed(futs), start=1):
                f = futs[fut]
                if stop.is_set():
                    for x in futs:
                        x.cancel()
                try:
                    for c, v in fut.result().items():
                        candidates[c].append(v)
                        tagged.append(pd.DataFrame({"channel": c, "well": f.well, "treatment": f.treatment or "",
                                                    "cell_type": f.cell_type or "", "b": v}))
                except Exception:  # noqa: BLE001 - the field's error is reported by the measure step
                    pass
                prog("measure", i, len(ready), "pooling spot brightness")
        if stop.is_set():
            raise InterruptedError("Stopped - nuclear masks are cached; run again to finish.")
    candidates = {c: (np.concatenate(v) if v else np.zeros(0)) for c, v in candidates.items()}
    spot_floors = resolve_spot_floors(params, speckle_ch, candidates) if do_spots else {}
    for c, d in spot_floors.items():
        sep = f", separability {d['separability']:.2f}" if np.isfinite(d["separability"]) else ""
        log(f"  {c}: spot brightness floor {d['floor']:.0f} ({d['mode']}{sep})")
        if np.isfinite(d["separability"]) and d["separability"] < SEPARABILITY_MIN:
            log(f"  WARNING {c}: spot brightness shows no clear dim/bright split - a fixed floor may remove real "
                "spots; check qc/spot_brightness.png or use the per-nucleus threshold.")
    floors = {c: d["floor"] for c, d in spot_floors.items()}
    # F-actin / phalloidin: readouts independent of the cell boundaries, with a bright cut-off from the reference
    from .actin import actin_channels

    actin_ch = actin_channels(ready, spot_channels, params)
    actin_cut, actin_info = {}, pd.DataFrame()
    if actin_ch:
        log(f"F-actin / phalloidin in {', '.join(actin_ch)}: whole-field and bright-structure readouts.")
        prog("measure", 0, len(ready), "F-actin bright cut-off")
        actin_cut, actin_info = resolve_actin_thresholds(ds, ready, actin_ch, params, masks, markers, log)
        actin_info.to_csv(adir / "actin_thresholds.csv", index=False)
    with ProcessPoolExecutor(max_workers=params.workers) as ex:
        futs = {ex.submit(_measure_task, f, _mask_path(masks, f), nuclear, spot_channels, params, ds.px_um,
                          ds.bit_ceiling, (qc_dir / f"{f.key}_qc.png") if params.qc_overlays else None,
                          floors, _cells_path(masks, f) if markers.get(f.key) else None, markers.get(f.key),
                          actin_cut.get(f.key)): f
                for f in ready}
        for i, fut in enumerate(as_completed(futs), start=1):
            f = futs[fut]
            if stop.is_set():
                for x in futs:
                    x.cancel()
            try:
                r = fut.result()
            except Exception as exc:  # noqa: BLE001 - e.g. a cancelled future
                r = FieldResult(f.key, error=str(exc) or "cancelled")
            if r.error:
                errors.append((f.key, r.error.splitlines()[0]))
                log(f"  {f.well} {f.field}: ERROR {r.error.splitlines()[0]}")
            else:
                results.append(r)
                n_ok = r.field_row.get("nuclei_qc_pass", 0)
                log(f"  {f.well} {f.field}: {len(r.nuclei)} nuclei ({n_ok} pass QC), {len(r.spots)} spots")
            prog("measure", i, len(ready), f"{f.well} {f.field}")
    if stop.is_set():
        raise InterruptedError("Stopped during measurement - nuclear masks are cached; run again to finish.")
    timings["measure_s"] = round(time.monotonic() - t_stage, 1)
    log(f"Measurement done in {timings['measure_s']} s.")
    t_stage = time.monotonic()

    # ---- stage 3: tables + statistics ----
    from .stats import (add_high_flags, condition_summary, default_orders, field_long_table, long_table,
                        replicate_summary, resolve_reference, run_tests, write_superplots)

    prog("stats", 0, 1, "tables")
    nuclei = pd.DataFrame([n for r in results for n in r.nuclei])
    for c in ("cell_type", "treatment", "antibody", "replicate"):
        if c in nuclei:
            nuclei[c] = nuclei[c].fillna("")
    # reference condition (e.g. WT · MOCK) and "high" cells relative to it
    reference, high_info = ("", ""), pd.DataFrame()
    if len(nuclei) and do_signal:
        o0 = default_orders(nuclei, params.treatment_order, params.cell_type_order)
        reference = resolve_reference(o0, params.reference_cell_type, params.reference_treatment)
        nuclei, high_info = add_high_flags(nuclei, signal_ch, reference, params.high_percentile)
        high_info.to_csv(adir / "high_thresholds.csv", index=False)
        log(f"Reference condition: {reference[0]} · {reference[1]} - 'high' cells are above its "
            f"{params.high_percentile:g}th percentile (high_thresholds.csv).")
    spots = pd.DataFrame([s for r in results for s in r.spots])
    coloc = pd.DataFrame([c for r in results for c in r.coloc])
    fields_df = pd.DataFrame([r.field_row for r in results])
    for df, name in ((nuclei, "nuclei"), (spots, "spots"), (coloc, "colocalization"), (fields_df, "fields")):
        df.to_csv(adir / f"{name}.csv", index=False)
    try:
        if not do_spots:
            raise RuntimeError("speckles not quantified")
        from .stats import write_spot_brightness_plot

        write_spot_brightness_plot(spots, pd.concat(tagged, ignore_index=True) if tagged else None, spot_floors,
                                   speckle_ch, ds.px_um, qc_dir / "spot_brightness.png")
        log("Wrote qc/spot_brightness.png (spot brightness per treatment, with any fixed floor).")
    except Exception as exc:  # noqa: BLE001 - a diagnostic plot never fails the run
        log(f"  spot brightness plot skipped: {exc}")
    log(f"Wrote nuclei.csv ({len(nuclei)}), spots.csv ({len(spots)}), colocalization.csv ({len(coloc)}).")

    for df in (nuclei, coloc):
        for c in ("cell_type", "treatment", "antibody", "replicate"):
            if c in df:
                df[c] = df[c].fillna("")
    long = long_table(nuclei, coloc, spot_channels, params.quantify_for) if len(nuclei) else pd.DataFrame()
    if len(long) and actin_ch and len(fields_df):  # per-field F-actin readouts (one value per field)
        long = pd.concat([long, field_long_table(fields_df, actin_ch)], ignore_index=True)
    rep = replicate_summary(long) if len(long) else pd.DataFrame()
    cond = condition_summary(rep) if len(rep) else pd.DataFrame()
    orders = default_orders(long, params.treatment_order, params.cell_type_order) if len(long) else {}
    if orders:
        log("Order: treatments " + ", ".join(orders["treatment"]) + "; cell types "
            + ", ".join(c or "-" for c in orders["cell_type"]))
    if orders:
        if not reference[0] and not reference[1]:
            reference = resolve_reference(orders, params.reference_cell_type, params.reference_treatment)
        orders["reference"] = reference
    omni, pairs = run_tests(rep, long, orders) if len(rep) else (pd.DataFrame(), pd.DataFrame())
    rep.to_csv(adir / "replicate_summary.csv", index=False)
    cond.to_csv(adir / "condition_summary.csv", index=False)
    omni.to_csv(adir / "stats_omnibus.csv", index=False)
    pairs.to_csv(adir / "stats_pairwise.csv", index=False)
    log("Wrote replicate_summary.csv, condition_summary.csv, stats_omnibus.csv, stats_pairwise.csv.")
    prog("stats", 0, 1, "SuperPlots")
    plots = write_superplots(long, rep, omni, adir / "superplots", orders, high_info) if len(rep) else []
    log(f"Wrote {len(plots)} summary figure(s) under superplots/ and spot_count_distribution.csv.")
    prog("stats", 1, 1, "done")
    timings["stats_s"] = round(time.monotonic() - t_stage, 1)

    manifest = {
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "software": _versions(),
        "input": str(Path(out_dir)),
        "nuclear_channel": nuclear, "spot_channels": spot_channels, "pixel_size_um": ds.px_um,
        "nucleus_method": method, "params": params.to_dict(), "spot_thresholds": spot_floors,
        "cell_method": cell_method, "cell_fields": len(cell_fields),
        "reference_condition": {"cell_type": reference[0], "treatment": reference[1],
                                "high_percentile": params.high_percentile},
        "cell_markers": sorted({markers[f.key] for f in cell_fields}),
        "actin_thresholds": actin_info.to_dict("records") if len(actin_info) else [],
        "fields_analysed": len(results), "errors": [{"field": k, "error": e} for k, e in errors],
        "nuclei_total": int(len(nuclei)), "nuclei_qc_pass": int(nuclei["qc_pass"].sum()) if len(nuclei) else 0,
        "replicate_unit": replicate_unit(params, ds),
        "seconds": round(time.monotonic() - t0, 1), "stage_seconds": timings,
    }
    (adir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    log(f"Done in {manifest['seconds']} s - {manifest['nuclei_qc_pass']} of {manifest['nuclei_total']} nuclei passed QC.")
    return AnalysisOutput(adir, nuclei, spots, coloc, cond, omni, pairs, plots, errors, method)


def preview_field(out_dir: "str | Path", key: str, params: AnalysisParams, png: Path) -> tuple[FieldResult, str]:
    """Segment + measure one field in this process (no cache) and write its QC overlay - for tuning."""
    ds = load_dataset(out_dir, params.replicate_column)
    fi = next(f for f in ds.fields if f.key == key)
    nuclear = params.nuclear_channel or ds.nuclear
    method = resolve_method(params)
    nuclei = segment_nuclei(tifffile.imread(fi.channels[nuclear]), params, ds.px_um, method)
    from .cells import cell_marker_for, segment_cells

    marker = cell_marker_for(fi.stains, fi.channels, params)
    cells = (segment_cells(tifffile.imread(fi.channels[marker]), tifffile.imread(fi.channels[nuclear]), nuclei, params,
                           ds.px_um) if marker else None)
    channels = params.spot_channels or ds.signal_channels
    floors = {c: float(v) for c, v in params.spot_fixed_threshold.items() if v and v > 0}
    auto = [c for c, v in params.spot_fixed_threshold.items() if not (v and v > 0)]
    man = Path(out_dir) / "analysis" / "analysis_manifest.json"
    last_run = json.loads(man.read_text(encoding="utf-8")) if man.exists() else {}
    if auto:  # an automatic floor needs every well: reuse the last full run's value
        last = last_run.get("spot_thresholds", {})
        floors.update({c: last[c]["floor"] for c in auto if c in last})
    from .actin import actin_channels

    actin = {}
    for c in actin_channels([fi], channels, params):  # the cut-off needs the reference wells: reuse the last run's
        if params.actin_bright_threshold and params.actin_bright_threshold > 0:
            actin[c] = float(params.actin_bright_threshold)
        else:
            hit = [x["threshold"] for x in last_run.get("actin_thresholds", [])
                   if x.get("channel") == c and x.get("antibody") == fi.antibody]
            actin[c] = float(hit[0]) if hit else float("nan")
    return measure_field(fi, nuclei, nuclear, channels, params, ds.px_um, ds.bit_ceiling, png, floors, cells,
                         marker, actin), method
