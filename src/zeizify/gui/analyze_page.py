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
        r = _row(c, "Flag nuclei by DAPI (× field median)", "dim / condensed (mitotic, apoptotic) nuclei are excluded from stats")
        self.qlo, self.qhi = _entry(r, d.qc_dapi_low, 55), _entry(r, d.qc_dapi_high, 55)
        ctk.CTkLabel(r, text="below", font=T.font(11), text_color=T.MUTED).pack(side="left", padx=(0, 4))
        self.qlo.pack(side="left")
        ctk.CTkLabel(r, text="or above", font=T.font(11), text_color=T.MUTED).pack(side="left", padx=6)
        self.qhi.pack(side="left")
        r = _row(c, "")
        self.border = _switch(r, "Exclude nuclei touching the image edge", d.exclude_border_nuclei)
        self.border.pack(side="left")
        self.gpu = _switch(r, "Use GPU (Cellpose)", d.use_gpu)
        self.gpu.pack(side="left", padx=20)
        ctk.CTkFrame(c, fg_color="transparent", height=8).pack()

        # spots
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Speckles / nuclear bodies",
                      "Detected inside each nucleus: white top-hat background removal, per-nucleus threshold with a "
                      "noise floor, size window").pack(fill="x", padx=18, pady=(14, 6))
        self.chan_row = _row(c, "Channels")
        self.chan_checks: dict[str, ctk.CTkCheckBox] = {}
        r = _row(c, "Spot diameter (µm)", "PML bodies ~0.3-1 µm · speckles ~0.5-3 µm")
        self.smin, self.smax = _entry(r, d.spot_min_diameter_um, 55), _entry(r, d.spot_max_diameter_um, 55)
        self.smin.pack(side="left")
        ctk.CTkLabel(r, text="to", font=T.font(11), text_color=T.MUTED).pack(side="left", padx=6)
        self.smax.pack(side="left")
        r = _row(c, "Noise floor (MADs)", "lower = more, fainter spots")
        self.sens_val = ctk.CTkLabel(r, text="", width=36, font=T.font(12, "bold"), text_color=T.TEXT)
        self.sens = ctk.CTkSlider(r, from_=2, to=12, number_of_steps=20, width=260, progress_color=T.ACCENT,
                                  button_color=T.ACCENT, button_hover_color=T.ACCENT_HOVER,
                                  command=lambda v: self.sens_val.configure(text=f"{v:.1f}"))
        self.sens.set(d.spot_sensitivity)
        self.sens_val.configure(text=f"{d.spot_sensitivity:.1f}")
        self.sens.pack(side="left")
        self.sens_val.pack(side="left", padx=8)
        r = _row(c, "")
        self.declump = _switch(r, "Split touching spots", d.spot_declump)
        self.declump.pack(side="left")
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
                                 fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER, border_color=T.BORDER)
            cb.select()
            cb.pack(side="left", padx=(0, 14))
            self.chan_checks[ch] = cb
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
                qc_dapi_low=float(self.qlo.get()), qc_dapi_high=float(self.qhi.get()),
                spot_channels=[c for c, cb in self.chan_checks.items() if cb.get()],
                spot_min_diameter_um=float(self.smin.get()), spot_max_diameter_um=float(self.smax.get()),
                spot_sensitivity=round(float(self.sens.get()), 1), spot_declump=bool(self.declump.get()),
                coloc_randomizations=int(self.rand.get()), replicate_column=REPLICATE_CHOICES[self.rep.get()],
                workers=max(1, int(self.workers.get())), qc_overlays=bool(self.qc.get()),
                treatment_order=_csv(self.t_order.get()), cell_type_order=_csv(self.ct_order.get()))
        except ValueError as exc:
            self.app.toast(f"Check the settings: {exc}", kind="warning")
            return None

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
        for pth in res.overview_plots[:4]:
            img = Image.open(pth)
            w = min(1000, img.width)
            ctk.CTkLabel(c, text="", image=ctk.CTkImage(img, img, size=(w, int(img.height * w / img.width)))).pack(
                padx=18, pady=(4, 14))

    def on_theme_change(self) -> None:
        pass
