"""
Analyze page - per-nucleus spot (speckle / nuclear body) and colocalization analysis of a
Zeizify output folder, with a one-field preview for tuning and a parallel full run.
"""

from __future__ import annotations

import threading
import time
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk
from PIL import Image

from . import theme as T
from .widgets import Card, ghost_button, primary_button, section_title

REPLICATE_CHOICES = {"Auto (plate-map Replicate column, else well)": "", "Well": "well",
                     "Field (technical - not independent)": "field"}


THRESHOLD_MODES = ["Per nucleus", "Fixed - automatic", "Fixed - value"]
CELL_MARKER_CHOICES = {"Auto (from the stain names)": ("auto", ""), "Off": ("off", "")}  # + one entry per channel
CELL_METHODS = {"Auto": "auto", "Cellpose-SAM": "cellpose", "Classical": "classical"}


def _csv(text: str) -> list[str]:
    return [t.strip() for t in text.split(",") if t.strip()]


def _entry(master, value, width=70) -> ctk.CTkEntry:
    e = ctk.CTkEntry(master, width=width, height=30, corner_radius=8, border_color=T.BORDER, fg_color=T.SURFACE_ALT,
                     font=T.font(12))
    e.insert(0, str(value))
    return e


def _row(master, label: str, hint: str = "") -> ctk.CTkFrame:
    r = ctk.CTkFrame(master, fg_color="transparent")
    r.pack(fill="x", padx=18, pady=3)
    ctk.CTkLabel(r, text=label, width=190, anchor="w", font=T.font(12), text_color=T.TEXT).pack(side="left")
    if hint:
        ctk.CTkLabel(r, text=hint, font=T.font(11), text_color=T.MUTED).pack(side="right")
    return r


def _switch(master, text, on=True) -> ctk.CTkSwitch:
    s = ctk.CTkSwitch(master, text=text, font=T.font(12), text_color=T.TEXT, progress_color=T.ACCENT)
    s.select() if on else s.deselect()
    return s


class AnalyzePage(ctk.CTkFrame):
    title = "Analyze"
    subtitle = "Nuclei, speckles / nuclear bodies and colocalization, per nucleus"

    def __init__(self, master, app):
        super().__init__(master, fg_color="transparent")
        from ..analysis.params import AnalysisParams
        from ..analysis.segment import cellpose_available

        self.app = app
        self.ds = None
        self.out: "Path | None" = None
        self.cancel = threading.Event()
        self.has_cellpose = cellpose_available()
        d = AnalysisParams()
        self.body = ctk.CTkScrollableFrame(self, fg_color="transparent", scrollbar_button_color=T.BORDER,
                                           scrollbar_button_hover_color=T.MUTED)
        self.body.pack(fill="both", expand=True)

        # input
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Input", "A Zeizify output folder - the one with intensity_report.csv").pack(fill="x", padx=18, pady=(14, 6))
        r = ctk.CTkFrame(c, fg_color="transparent")
        r.pack(fill="x", padx=18, pady=(0, 6))
        self.folder_btn = ctk.CTkButton(r, text="▢  Choose output folder…", anchor="w", height=32, corner_radius=16,
                                        fg_color=T.SURFACE_ALT, hover_color=T.BORDER, text_color=T.MUTED,
                                        font=T.font(12), command=self.pick_folder)
        self.folder_btn.pack(side="left", fill="x", expand=True)
        self.info = ctk.CTkLabel(c, text="", font=T.font(12), text_color=T.MUTED, anchor="w", justify="left", wraplength=900)
        self.info.pack(fill="x", padx=18, pady=(0, 14))

        # nuclei
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Nuclei", "Segmented from the nuclear stain").pack(fill="x", padx=18, pady=(14, 6))
        r = _row(c, "Method", "Cellpose: Stringer et al. 2021 · Classical: Otsu + watershed")
        values = ["Auto", "Cellpose", "Classical"]
        self.method = ctk.CTkSegmentedButton(r, values=values, font=T.font(12), selected_color=T.ACCENT,
                                             selected_hover_color=T.ACCENT_HOVER, fg_color=T.SURFACE_ALT,
                                             unselected_color=T.SURFACE_ALT, unselected_hover_color=T.BORDER,
                                             text_color=T.TEXT)
        self.method.set("Auto")
        self.method.pack(side="left")
        if not self.has_cellpose:
            ctk.CTkLabel(c, text="Cellpose is not installed - Auto uses the classical method "
                                 "(pip install cellpose; see README for GPU setup).",
                         font=T.font(11), text_color=T.WARNING, anchor="w").pack(fill="x", padx=18)
        r = _row(c, "Typical nucleus diameter (µm)")
        self.diam = _entry(r, d.nucleus_diameter_um)
        self.diam.pack(side="left")
        r = _row(c, "")
        self.border = _switch(r, "Exclude nuclei touching the image edge", d.exclude_border_nuclei)
        self.border.pack(side="left")
        self.gpu = _switch(r, "Use GPU (Cellpose)", d.use_gpu)
        self.gpu.pack(side="left", padx=20)
        ctk.CTkFrame(c, fg_color="transparent", height=8).pack()

        # cells / cytoplasm
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Cells and cytoplasm",
                      "Whole-cell boundaries from a cell-body / edge stain such as F-actin, for wells that have one. "
                      "Cytoplasm = cell minus nucleus; adds cytoplasmic intensity and nuclear / cytoplasmic ratios."
                      ).pack(fill="x", padx=18, pady=(14, 6))
        r = _row(c, "Boundary stain", "Auto: the channel whose stain is F-actin, phalloidin, CellMask, tubulin or WGA")
        self.cell_marker = ctk.CTkOptionMenu(r, values=list(CELL_MARKER_CHOICES), width=260, height=30, corner_radius=8,
                                             font=T.font(12), fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT,
                                             button_hover_color=T.BORDER, text_color=T.TEXT, dropdown_font=T.font(12))
        self.cell_marker.set(next(iter(CELL_MARKER_CHOICES)))
        self.cell_marker.pack(side="left")
        r = _row(c, "Method", "Cellpose-SAM separates touching cells; Classical grows each nucleus outwards")
        self.cell_method = ctk.CTkSegmentedButton(r, values=list(CELL_METHODS), font=T.font(12),
                                                  selected_color=T.ACCENT, selected_hover_color=T.ACCENT_HOVER,
                                                  fg_color=T.SURFACE_ALT, unselected_color=T.SURFACE_ALT,
                                                  unselected_hover_color=T.BORDER, text_color=T.TEXT)
        self.cell_method.set(next(iter(CELL_METHODS)))
        self.cell_method.pack(side="left")
        r = _row(c, "Typical cell diameter (µm)", "sets the scale for both methods")
        self.cell_diam = _entry(r, d.cell_diameter_um)
        self.cell_diam.pack(side="left")
        r = _row(c, "")
        self.cell_border = _switch(r, "No cytoplasm for cells cut by the image edge", d.exclude_border_cells)
        self.cell_border.pack(side="left")
        self.cell_info = ctk.CTkLabel(c, text="", font=T.font(11), text_color=T.MUTED, anchor="w", justify="left",
                                      wraplength=900)
        self.cell_info.pack(fill="x", padx=18, pady=(2, 12))

        # spots
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Speckles / nuclear bodies",
                      "Detected inside each nucleus: white top-hat background removal, per-nucleus threshold with a "
                      "noise floor, size window").pack(fill="x", padx=18, pady=(14, 6))
        self.chan_row = _row(c, "Channels")
        ctk.CTkLabel(c, text="QUANTIFY PER FLUOROPHORE - Speckles: spot detection, counts, sizes, overlap · Overall "
                             "signal: mean and total intensity of nucleus / cytoplasm / cell, speckle signal included",
                     font=T.font(10, "bold"), text_color=T.ACCENT, anchor="w", justify="left", wraplength=900).pack(
            fill="x", padx=18, pady=(8, 0))
        self.quant_box = ctk.CTkFrame(c, fg_color="transparent")
        self.quant_box.pack(fill="x")
        self.quant_checks: dict[str, tuple] = {}
        self.quant_hint = ctk.CTkLabel(self.quant_box, text="Choose an output folder to list its fluorophores.",
                                       font=T.font(11), text_color=T.MUTED, anchor="w")
        self.quant_hint.pack(fill="x", padx=18)
        self.chan_checks: dict[str, ctk.CTkCheckBox] = {}
        self.thr_box = ctk.CTkFrame(c, fg_color="transparent")
        self.thr_box.pack(fill="x")
        self.thr_widgets: dict[str, tuple] = {}
        self.spot_note = ctk.CTkLabel(c, text="Spot options are off: no fluorophore has Speckles ticked.",
                                      font=T.font(11), text_color=T.MUTED, anchor="w")
        r = _row(c, "Spot diameter (µm)", "PML bodies ~0.3-1 µm · speckles ~0.5-3 µm")
        self.spot_rows = [r]
        self.smin, self.smax = _entry(r, d.spot_min_diameter_um, 55), _entry(r, d.spot_max_diameter_um, 55)
        self.smin.pack(side="left")
        ctk.CTkLabel(r, text="to", font=T.font(11), text_color=T.MUTED).pack(side="left", padx=6)
        self.smax.pack(side="left")
        r = _row(c, "Noise floor (MADs)", "lower = more, fainter spots")
        self.spot_rows.append(r)
        self.sens_val = ctk.CTkLabel(r, text="", width=36, font=T.font(12, "bold"), text_color=T.TEXT)
        self.sens = ctk.CTkSlider(r, from_=2, to=12, number_of_steps=20, width=260, progress_color=T.ACCENT,
                                  button_color=T.ACCENT, button_hover_color=T.ACCENT_HOVER,
                                  command=lambda v: self.sens_val.configure(text=f"{v:.1f}"))
        self.sens.set(d.spot_sensitivity)
        self.sens_val.configure(text=f"{d.spot_sensitivity:.1f}")
        self.sens.pack(side="left")
        self.sens_val.pack(side="left", padx=8)
        r = _row(c, "")
        self.spot_rows.append(r)
        self.declump = _switch(r, "Split touching spots", d.spot_declump)
        self.declump.pack(side="left")
        r = _row(c, "Quantify spots in", "tick both for nucleus, cytoplasm and whole-cell graphs; cytoplasm needs "
                                         "cell boundaries (e.g. F-actin)")
        self.spot_rows.append(r)
        self.in_nucleus = ctk.CTkCheckBox(r, text="Nucleus", font=T.font(12), text_color=T.TEXT, checkbox_width=18,
                                          checkbox_height=18, corner_radius=5, fg_color=T.ACCENT,
                                          hover_color=T.ACCENT_HOVER, border_color=T.BORDER)
        self.in_cyto = ctk.CTkCheckBox(r, text="Cytoplasm", font=T.font(12), text_color=T.TEXT, checkbox_width=18,
                                       checkbox_height=18, corner_radius=5, fg_color=T.ACCENT,
                                       hover_color=T.ACCENT_HOVER, border_color=T.BORDER)
        self.in_nucleus.select()
        self.in_nucleus.pack(side="left", padx=(0, 14))
        self.in_cyto.pack(side="left")
        ctk.CTkFrame(c, fg_color="transparent", height=8).pack()

        # coloc + run settings
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Colocalization & statistics",
                      "Every pair of selected channels, per nucleus: Pearson, Spearman, Li ICQ, Manders tM1/tM2 (Costes), "
                      "and spot overlap. Tests use replicate medians (SuperPlots).").pack(fill="x", padx=18, pady=(14, 6))
        r = _row(c, "Costes randomizations", "0 = off; 100 adds a significance test per nucleus (slower)")
        self.rand = _entry(r, d.coloc_randomizations)
        self.rand.pack(side="left")
        r = _row(c, "Replicate unit")
        self.rep = ctk.CTkOptionMenu(r, values=list(REPLICATE_CHOICES), width=320, height=30, corner_radius=8,
                                     font=T.font(12), fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT,
                                     button_hover_color=T.BORDER, text_color=T.TEXT, dropdown_font=T.font(12))
        self.rep.pack(side="left")
        r = _row(c, "Treatment order", "comma-separated; blank = control first, then plate-map order")
        self.t_order = _entry(r, "", 420)
        self.t_order.pack(side="left")
        r = _row(c, "Cell-type order", "comma-separated; blank = wild type first - all cell types share each graph")
        self.ct_order = _entry(r, "", 420)
        self.ct_order.pack(side="left")
        r = _row(c, "Reference condition", "every condition is compared with it; it also sets the 'high' cut-off")
        self.ref_ct = ctk.CTkOptionMenu(r, values=["Auto (wild type)"], width=180, height=30, corner_radius=8,
                                        font=T.font(12), fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT,
                                        button_hover_color=T.BORDER, text_color=T.TEXT, dropdown_font=T.font(12))
        self.ref_ct.set("Auto (wild type)")
        self.ref_ct.pack(side="left")
        self.ref_t = ctk.CTkOptionMenu(r, values=["Auto (control)"], width=180, height=30, corner_radius=8,
                                       font=T.font(12), fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT,
                                       button_hover_color=T.BORDER, text_color=T.TEXT, dropdown_font=T.font(12))
        self.ref_t.set("Auto (control)")
        self.ref_t.pack(side="left", padx=8)
        r = _row(c, "'High' cells", "above this percentile of the reference condition's cells")
        self.high_pct = _entry(r, d.high_percentile, 60)
        self.high_pct.pack(side="left")
        ctk.CTkLabel(r, text="th percentile", font=T.font(11), text_color=T.MUTED).pack(side="left", padx=6)
        r = _row(c, "Parallel workers", "CPU processes for measurement")
        self.workers = _entry(r, d.workers)
        self.workers.pack(side="left")
        r = _row(c, "")
        self.qc = _switch(r, "Save QC overlay per field", d.qc_overlays)
        self.qc.pack(side="left")
        ctk.CTkFrame(c, fg_color="transparent", height=8).pack()

        # preview
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Preview", "Run the current settings on one field to check segmentation before the full run. "
                                    "White = analysed nucleus, red = QC-flagged, green/magenta = spots.").pack(
            fill="x", padx=18, pady=(14, 6))
        r = ctk.CTkFrame(c, fg_color="transparent")
        r.pack(fill="x", padx=18, pady=(0, 6))
        self.field_menu = ctk.CTkOptionMenu(r, values=["(load a folder)"], width=420, height=30, corner_radius=8,
                                            font=T.font(12), fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT,
                                            button_hover_color=T.BORDER, text_color=T.TEXT, dropdown_font=T.font(12))
        self.field_menu.pack(side="left")
        self.preview_btn = ghost_button(r, "Preview", self.preview, width=90)
        self.preview_btn.pack(side="left", padx=8)
        self.preview_info = ctk.CTkLabel(c, text="", font=T.font(12), text_color=T.MUTED, anchor="w", justify="left",
                                         wraplength=900)
        self.preview_info.pack(fill="x", padx=18)
        self.preview_img = ctk.CTkLabel(c, text="")
        self.preview_img.pack(padx=18, pady=(4, 14))

        # run
        self.progress_slot = ctk.CTkFrame(self.body, fg_color="transparent", height=1)
        self.progress_slot.pack(fill="x")
        bar = Card(self.body)
        bar.pack(fill="x", pady=(0, 8))
        self.status = ctk.CTkLabel(bar, text="Results are written to <output folder>/analysis/.", font=T.font(11),
                                   text_color=T.MUTED, anchor="w", justify="left", wraplength=640)
        self.status.pack(side="left", padx=18, pady=14, fill="x", expand=True)
        self.run_btn = primary_button(bar, "Run analysis", self.run, width=140)
        self.run_btn.pack(side="right", padx=14, pady=12)
        self.results_slot = ctk.CTkFrame(self.body, fg_color="transparent", height=1)
        self.results_slot.pack(fill="x")

    # ---------- input ----------

    def pick_folder(self) -> None:
        p = filedialog.askdirectory(title="Select a Zeizify output folder (with intensity_report.csv)")
        if p:
            self.load(Path(p))

    def load(self, out: Path) -> None:
        from ..analysis.io import load_dataset

        try:
            ds = load_dataset(out)
        except Exception as exc:  # noqa: BLE001
            self.info.configure(text=f"Could not load: {exc}", text_color=T.DANGER)
            return
        self.ds, self.out = ds, out
        self.folder_btn.configure(text=f"▢  {out}", text_color=T.TEXT)
        wells = sorted({f.well for f in ds.fields}, key=lambda w: (w[0], int(w[1:])))
        chans = ", ".join(f"{c} ({', '.join(ds.stains_for(c))})" if ds.stains_for(c) else c for c in ds.channels)
        px = f"{ds.px_um:.4f} µm/px" if ds.px_um else "pixel size unknown"
        self.info.configure(text=f"{len(ds.fields)} fields · {len(wells)} wells · channels: {chans} · nuclear: "
                                 f"{ds.nuclear or '?'} · {px}", text_color=T.MUTED)
        for w in self.chan_row.winfo_children()[1:]:
            w.destroy()
        self.chan_checks = {}
        for ch in ds.signal_channels:
            st = ", ".join(ds.stains_for(ch))
            cb = ctk.CTkCheckBox(self.chan_row, text=f"{ch}" + (f" · {st}" if st else ""), font=T.font(12),
                                 text_color=T.TEXT, checkbox_width=18, checkbox_height=18, corner_radius=5,
                                 fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER, border_color=T.BORDER,
                                 command=self._sync_spot_controls)
            cb.select()
            cb.pack(side="left", padx=(0, 14))
            self.chan_checks[ch] = cb
        for w in self.thr_box.winfo_children():
            w.destroy()
        self.thr_widgets = {}
        for ch in ds.signal_channels:
            st = ", ".join(ds.stains_for(ch))
            r = _row(self.thr_box, f"Threshold · {ch}" + (f" ({st})" if st else ""),
                     "Fixed: one brightness floor for every well - for stains whose bodies vanish in some "
                     "conditions" if ch == ds.signal_channels[0] else "")
            menu = ctk.CTkOptionMenu(r, values=THRESHOLD_MODES, width=200, height=30, corner_radius=8, font=T.font(12),
                                     fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT, button_hover_color=T.BORDER,
                                     text_color=T.TEXT, dropdown_font=T.font(12))
            menu.set(THRESHOLD_MODES[0])
            menu.pack(side="left")
            val = _entry(r, "", 90)
            val.pack(side="left", padx=8)
            ctk.CTkLabel(r, text="value (Fixed - value only)", font=T.font(11), text_color=T.MUTED).pack(side="left")
            self.thr_widgets[ch] = (menu, val, r)
        values = list(CELL_MARKER_CHOICES) + [f"{ch} · {', '.join(ds.stains_for(ch))}" if ds.stains_for(ch) else ch
                                              for ch in ds.signal_channels]
        self.cell_marker.configure(values=values)
        if self.cell_marker.get() not in values:
            self.cell_marker.set(values[0])
        self.cell_info.configure(text=self._describe_cells(ds))
        from ..analysis.stats import is_control, is_wildtype, level_order
        import pandas as pd

        meta = pd.DataFrame([{"cell_type": f.cell_type, "treatment": f.treatment, "well": f.well} for f in ds.fields])
        cts = level_order(meta, "cell_type", first=is_wildtype)
        trs = level_order(meta, "treatment", first=is_control)
        for menu, auto, values in ((self.ref_ct, f"Auto (wild type: {cts[0] if cts else '-'})", cts),
                                   (self.ref_t, f"Auto (control: {trs[0] if trs else '-'})", trs)):
            menu.configure(values=[auto, *values])
            if menu.get() not in values:
                menu.set(auto)
        self._build_quant_rows(ds)
        self.field_menu.configure(values=[f.key for f in ds.fields])
        self.field_menu.set(ds.fields[0].key if ds.fields else "")

    def on_show(self) -> None:
        run = self.app.pages["run"]
        if self.ds is None and getattr(run, "out", None) and (Path(run.out) / "intensity_report.csv").exists():
            self.load(Path(run.out))

    # ---------- params ----------

    def params(self):
        from ..analysis.params import AnalysisParams

        try:
            return AnalysisParams(
                nucleus_method=self.method.get().lower(), nucleus_diameter_um=float(self.diam.get()),
                exclude_border_nuclei=bool(self.border.get()), use_gpu=bool(self.gpu.get()),
                spot_channels=[c for c, cb in self.chan_checks.items() if cb.get()],
                spot_min_diameter_um=float(self.smin.get()), spot_max_diameter_um=float(self.smax.get()),
                spot_sensitivity=round(float(self.sens.get()), 1), spot_declump=bool(self.declump.get()),
                spot_compartments=self._compartments(), quantify=self._quantify()[0],
                channel_quantify=self._quantify()[1],
                coloc_randomizations=int(self.rand.get()), replicate_column=REPLICATE_CHOICES[self.rep.get()],
                workers=max(1, int(self.workers.get())), qc_overlays=bool(self.qc.get()),
                cell_segmentation=self._marker_choice()[0], cell_marker=self._marker_choice()[1],
                reference_cell_type="" if self.ref_ct.get().startswith("Auto") else self.ref_ct.get(),
                reference_treatment="" if self.ref_t.get().startswith("Auto") else self.ref_t.get(),
                high_percentile=float(self.high_pct.get()),
                cell_method=CELL_METHODS[self.cell_method.get()], cell_diameter_um=float(self.cell_diam.get()),
                exclude_border_cells=bool(self.cell_border.get()),
                treatment_order=_csv(self.t_order.get()), cell_type_order=_csv(self.ct_order.get()),
                spot_fixed_threshold=self._fixed_thresholds())
        except ValueError as exc:
            self.app.toast(f"Check the settings: {exc}", kind="warning")
            return None

    def _quantify(self) -> tuple[list[str], dict[str, list[str]]]:
        """(everything quantified in any channel, per-channel choice) for the analysed channels."""
        per = {}
        for ch, (spk, sig) in self.quant_checks.items():
            if ch in self.chan_checks and not self.chan_checks[ch].get():
                continue  # channel not analysed
            q = [name for name, cb in (("speckles", spk), ("signal", sig)) if cb.get()]
            if not q:
                raise ValueError(f"{ch}: tick Speckles and/or Overall signal (or untick the channel)")
            per[ch] = q
        union = sorted({x for q in per.values() for x in q}, key=["speckles", "signal"].index) or ["speckles", "signal"]
        return union, per

    def _compartments(self) -> list[str]:
        if "speckles" not in self._quantify()[0]:
            return ["nucleus"]  # no spot detection anywhere: the (disabled) choice does not matter
        out = [name for name, cb in (("nucleus", self.in_nucleus), ("cytoplasm", self.in_cyto)) if cb.get()]
        if not out:
            raise ValueError("tick Nucleus and/or Cytoplasm under 'Quantify spots in'")
        if "cytoplasm" in out and self._marker_choice()[0] == "off":
            raise ValueError("cytoplasmic spots need cell boundaries - set Boundary stain to Auto or a channel")
        return out

    def _build_quant_rows(self, ds) -> None:
        """One row per fluorophore: Speckles / Overall signal. A cell-boundary stain (F-actin, ...) starts with
        Overall signal only - it has no speckles to count."""
        from ..analysis.cells import CELL_MARKER_KEYWORDS

        # keep the choices when the same folder is reloaded; a different folder starts from the defaults
        same = getattr(self, "_quant_dir", None) == ds.out_dir
        old = {ch: (spk.get(), sig.get()) for ch, (spk, sig) in self.quant_checks.items()} if same else {}
        self._quant_dir = ds.out_dir
        for w in self.quant_box.winfo_children():
            w.destroy()
        self.quant_checks = {}
        for ch in ds.signal_channels:
            stains = ", ".join(ds.stains_for(ch))
            r = _row(self.quant_box, f"{ch}" + (f" · {stains}" if stains else ""))
            spk = ctk.CTkCheckBox(r, text="Speckles", font=T.font(12), text_color=T.TEXT, checkbox_width=18,
                                  checkbox_height=18, corner_radius=5, fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER,
                                  border_color=T.BORDER, command=self._sync_spot_controls)
            sig = ctk.CTkCheckBox(r, text="Overall signal", font=T.font(12), text_color=T.TEXT, checkbox_width=18,
                                  checkbox_height=18, corner_radius=5, fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER,
                                  border_color=T.BORDER)
            marker = any(k in stains.lower() for k in CELL_MARKER_KEYWORDS)
            on_spk, on_sig = old.get(ch, (not marker, True))
            spk.select() if on_spk else spk.deselect()
            sig.select() if on_sig else sig.deselect()
            spk.pack(side="left", padx=(0, 14))
            sig.pack(side="left")
            if marker and ch not in old:
                ctk.CTkLabel(r, text="cell-boundary stain: overall signal only by default", font=T.font(11),
                             text_color=T.MUTED).pack(side="left", padx=12)
            self.quant_checks[ch] = (spk, sig)
        self._sync_spot_controls()

    def _speckle_channels(self) -> set[str]:
        """Analysed channels with Speckles ticked."""
        return {ch for ch, (spk, _sig) in self.quant_checks.items()
                if spk.get() and (ch not in self.chan_checks or self.chan_checks[ch].get())}

    @staticmethod
    def _enable_row(row, on: bool) -> None:
        """Enable / grey out every control in a settings row; the row label is dimmed with it."""
        kids = row.winfo_children()
        for i, w in enumerate(kids):
            if isinstance(w, ctk.CTkLabel):
                if i == 0:  # the row label (hints keep their muted colour)
                    w.configure(text_color=T.TEXT if on else T.MUTED)
                continue
            try:
                w.configure(state="normal" if on else "disabled")
            except (ValueError, TypeError, AttributeError):
                pass

    def _sync_spot_controls(self) -> None:
        """Spot options apply only where speckles are detected: a fluorophore's threshold row needs its own
        Speckles box; diameter, noise floor, splitting and 'Quantify spots in' need Speckles in at least one."""
        speckle = self._speckle_channels()
        for ch, (spk, sig) in self.quant_checks.items():
            analysed = ch not in self.chan_checks or bool(self.chan_checks[ch].get())
            for cb in (spk, sig):
                cb.configure(state="normal" if analysed else "disabled")
        for ch, (_menu, _val, row) in self.thr_widgets.items():
            self._enable_row(row, ch in speckle)
        for row in self.spot_rows:
            self._enable_row(row, bool(speckle))
        if speckle or not self.quant_checks:
            self.spot_note.pack_forget()
        else:
            self.spot_note.pack(fill="x", padx=18, pady=(4, 0), before=self.spot_rows[0])

    def _marker_choice(self) -> tuple[str, str]:
        v = self.cell_marker.get()
        if v in CELL_MARKER_CHOICES:
            return CELL_MARKER_CHOICES[v]
        return "auto", v.split(" ")[0]  # "AF488 · F-Actin" -> a fixed channel

    def _describe_cells(self, ds) -> str:
        """Which staining panels get cell boundaries, from which channel (default settings)."""
        from ..analysis.cells import cell_marker_for
        from ..analysis.params import AnalysisParams

        p = AnalysisParams()
        by_panel = {}
        for f in ds.fields:
            m = cell_marker_for(f.stains, f.channels, p)
            by_panel.setdefault(f.antibody or "(no panel)", set()).add(f"{m} ({f.stains.get(m, '')})" if m else "")
        on = [f"{ab}: {', '.join(sorted(s - {''}))}" for ab, s in by_panel.items() if s - {""}]
        off = [ab for ab, s in by_panel.items() if not (s - {""})]
        if not on:
            return "Auto: no well has a cell-boundary stain - only nuclei are measured. Choose a channel above to force it."
        return ("Auto finds cell boundaries in " + "; ".join(on)
                + (f".  Nuclei only: {', '.join(off)}." if off else "."))

    def _fixed_thresholds(self) -> dict[str, float]:
        out = {}
        speckle = self._speckle_channels()
        for ch, (menu, val, _r) in self.thr_widgets.items():
            if ch not in speckle:
                continue  # no spot detection in this channel: its threshold is not used
            mode = menu.get()
            if mode == "Fixed - automatic":
                out[ch] = 0.0
            elif mode == "Fixed - value":
                v = float(val.get())  # ValueError -> "Check the settings"
                if v <= 0:
                    raise ValueError(f"{ch}: the fixed threshold must be > 0")
                out[ch] = v
        return out

    def _ready(self) -> bool:
        if self.ds is None:
            self.app.toast("Choose a Zeizify output folder first.", kind="warning")
            return False
        if not any(cb.get() for cb in self.chan_checks.values()):
            self.app.toast("Select at least one spot channel.", kind="warning")
            return False
        return True

    # ---------- preview ----------

    def preview(self) -> None:
        if not self._ready():
            return
        p = self.params()
        if p is None:
            return
        key = self.field_menu.get()
        png = self.out / "analysis" / "preview" / f"{key}_preview.png"
        self.preview_btn.configure(state="disabled", text="Working…")
        self.preview_info.configure(text="Segmenting and measuring one field…", text_color=T.MUTED)
        t0 = time.monotonic()

        def work():
            from ..analysis.run import preview_field

            return preview_field(self.out, key, p, png)

        def done(res):
            import numpy as np

            r, method = res
            self.preview_btn.configure(state="normal", text="Preview")
            ok = [n for n in r.nuclei if n["qc_pass"]]
            parts = [f"{len(r.nuclei)} nuclei ({len(ok)} pass QC) · {method} · {time.monotonic() - t0:.1f} s"]
            for c in p.spot_channels:
                v = [n[f"{c}_spots_n"] for n in ok]
                if v:
                    parts.append(f"{c}: median {np.median(v):.0f} spots/nucleus (range {min(v)}-{max(v)})")
            cl = [x["pearson"] for x in r.coloc if x["qc_pass"]]
            if cl:
                parts.append(f"median Pearson {np.nanmedian(cl):.2f}")
            reasons = sorted({n["qc_reason"] for n in r.nuclei if n["qc_reason"]})
            if reasons:
                parts.append("flagged: " + "; ".join(reasons))
            self.preview_info.configure(text="   ·   ".join(parts), text_color=T.TEXT)
            img = Image.open(png)
            w = min(900, img.width)
            self.preview_img.configure(image=ctk.CTkImage(img, img, size=(w, int(img.height * w / img.width))))

        def fail(exc):
            self.preview_btn.configure(state="normal", text="Preview")
            self.preview_info.configure(text=f"Preview failed: {exc}", text_color=T.DANGER)

        self.app.run_bg(work, done, fail)

    # ---------- full run ----------

    def run(self) -> None:
        if not self._ready():
            return
        p = self.params()
        if p is None:
            return
        from .contrast_page import ExportPanel

        self.cancel = threading.Event()
        for w in self.progress_slot.winfo_children() + self.results_slot.winfo_children():
            w.destroy()
        panel = ExportPanel(self.progress_slot, on_stop=lambda: (
            self.cancel.set(), panel.stop_btn.configure(text="Stopping…", state="disabled")))
        panel.STEPS = ["Nuclei", "Measure", "Statistics"]
        panel.title.configure(text="Analysing…")
        panel.counts.pack_forget()
        panel.pack(fill="x", pady=(0, 12))
        panel.set_step(0)
        self.run_btn.configure(state="disabled", text="Running…", fg_color=T.BORDER, text_color=T.MUTED)
        self.status.configure(text="Running - progress is shown above. You can keep using the rest of the app.",
                              text_color=T.MUTED)
        self.after(50, lambda: self.body._parent_canvas.yview_moveto(1.0))
        stage_ix = {"nuclei": 0, "measure": 1, "stats": 2}
        current = {"stage": "nuclei"}

        def progress(stage, done, total, item):
            def ui():
                if stage != current["stage"]:
                    current["stage"] = stage
                    panel.set_step(stage_ix[stage])
                if stage == "stats":
                    panel.current.configure(text=f"{item}")
                else:
                    panel.progress(done, total, item)
            self.app.post(ui)

        def work():
            from ..analysis.run import run_analysis

            return run_analysis(self.out, p, log=lambda s: self.app.post(panel.add_log, s), progress=progress,
                                cancel=self.cancel)

        def finish():
            self.run_btn.configure(state="normal", text="Run analysis", fg_color=T.ACCENT, text_color=T.ON_ACCENT)

        def done(res):
            finish()
            panel.finish("Analysis complete", "success")
            n_ok = int(res.nuclei["qc_pass"].sum()) if len(res.nuclei) else 0
            self.status.configure(text=f"{len(res.nuclei)} nuclei ({n_ok} passed QC) · results in {res.analysis_dir}",
                                  text_color=T.SUCCESS)
            self._results(res)
            self.app.toast("Analysis complete.", kind="success")

        def fail(exc):
            finish()
            stopped = isinstance(exc, InterruptedError)
            panel.finish("Analysis stopped" if stopped else "Analysis failed", "warning" if stopped else "error")
            panel.add_log(("" if stopped else "ERROR: ") + str(exc))
            self.status.configure(text=str(exc), text_color=T.WARNING if stopped else T.DANGER)

        self.app.run_bg(work, done, fail)

    def _results(self, res) -> None:
        from .run_page import open_path

        c = Card(self.results_slot)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Results", "nuclei.csv · spots.csv · colocalization.csv · condition_summary.csv · "
                                    "stats_omnibus.csv · stats_pairwise.csv · superplots/ · qc/").pack(
            fill="x", padx=18, pady=(14, 6))
        if len(res.omnibus) and res.omnibus["welch_p"].isna().all():
            ctk.CTkLabel(c, text="No p-values: each treatment needs at least 2 replicates (wells, or a Replicate column "
                                 "in the plate map). Effect sizes are still reported in stats_pairwise.csv.",
                         font=T.font(12), text_color=T.WARNING, anchor="w", justify="left", wraplength=900).pack(
                fill="x", padx=18)
        if res.errors:
            ctk.CTkLabel(c, text=f"{len(res.errors)} field(s) failed: " + "; ".join(k for k, _ in res.errors[:5]),
                         font=T.font(12), text_color=T.DANGER, anchor="w").pack(fill="x", padx=18)
        a = ctk.CTkFrame(c, fg_color="transparent")
        a.pack(fill="x", padx=18, pady=8)
        primary_button(a, "Open results folder", lambda: open_path(res.analysis_dir), width=160).pack(side="left")
        ghost_button(a, "Open SuperPlots", lambda: open_path(res.analysis_dir / "superplots"), width=130).pack(
            side="left", padx=6)
        ghost_button(a, "Open QC overlays", lambda: open_path(res.analysis_dir / "qc"), width=130).pack(side="left")
        ghost_button(a, "Single-cell figure  ▣", self._to_figures, width=160).pack(side="left", padx=6)
        for pth in res.overview_plots[:4]:
            img = Image.open(pth)
            w = min(1000, img.width)
            ctk.CTkLabel(c, text="", image=ctk.CTkImage(img, img, size=(w, int(img.height * w / img.width)))).pack(
                padx=18, pady=(4, 14))

    def _to_figures(self) -> None:
        page = self.app.pages["figure"]
        if self.out is not None:
            page.load(self.out)
        self.app.show("figure")

    def on_theme_change(self) -> None:
        pass
