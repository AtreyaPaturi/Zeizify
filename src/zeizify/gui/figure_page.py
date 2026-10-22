"""
Figures page - publication single-cell figures: one representative nucleus per cell type x
treatment (picked automatically from the Analyze results, or chosen from ranked candidates),
cropped to a custom size and laid out with the stain panels as rows or columns.
"""

from __future__ import annotations

import io
import threading
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk
import numpy as np
from PIL import Image

from . import theme as T
from .widgets import Card, ghost_button, primary_button, section_title

ROWS = {"Treatments as rows": "treatment", "Stains as rows": "panel"}
BLOCKS = {"Cell types side by side": "side", "Cell types stacked": "stacked"}
CONTRAST = {"Same range in every condition": "shared", "Each panel separately": "panel"}
OUTLINE = {"Nuclear-stain panels": "nuclear", "All panels": "all", "None": "none"}
SCALEBAR = {"Bottom-right panel": "last", "Every panel": "all", "None": "none"}
FONTS = ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"]
CUSTOM_WIDTH = "Custom width (mm)"


def _entry(master, value="", width=70) -> ctk.CTkEntry:
    e = ctk.CTkEntry(master, width=width, height=30, corner_radius=8, border_color=T.BORDER, fg_color=T.SURFACE_ALT,
                     font=T.font(12))
    if value != "":
        e.insert(0, str(value))
    return e


def _row(master, label: str, hint: str = "") -> ctk.CTkFrame:
    r = ctk.CTkFrame(master, fg_color="transparent")
    r.pack(fill="x", padx=18, pady=3)
    ctk.CTkLabel(r, text=label, width=170, anchor="w", font=T.font(12), text_color=T.TEXT).pack(side="left")
    if hint:
        ctk.CTkLabel(r, text=hint, font=T.font(11), text_color=T.MUTED).pack(side="right")
    return r


def _menu(master, values, width=230, command=None) -> ctk.CTkOptionMenu:
    m = ctk.CTkOptionMenu(master, values=list(values), width=width, height=30, corner_radius=8, font=T.font(12),
                          fg_color=T.SURFACE_ALT, button_color=T.SURFACE_ALT, button_hover_color=T.BORDER,
                          text_color=T.TEXT, dropdown_font=T.font(12), command=command)
    m.set(list(values)[0])
    return m


def _seg(master, values) -> ctk.CTkSegmentedButton:
    s = ctk.CTkSegmentedButton(master, values=list(values), font=T.font(12), selected_color=T.ACCENT,
                               selected_hover_color=T.ACCENT_HOVER, fg_color=T.SURFACE_ALT,
                               unselected_color=T.SURFACE_ALT, unselected_hover_color=T.BORDER, text_color=T.TEXT)
    s.set(list(values)[0])
    return s


def _check(master, text, on=True, command=None) -> ctk.CTkCheckBox:
    cb = ctk.CTkCheckBox(master, text=text, font=T.font(12), text_color=T.TEXT, checkbox_width=18, checkbox_height=18,
                         corner_radius=5, fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER, border_color=T.BORDER,
                         command=command)
    cb.select() if on else cb.deselect()
    return cb


class FieldPicker(ctk.CTkToplevel):
    """Whole-field viewer for one condition: click a nucleus to choose it for the figure."""

    MAX_W, MAX_H = 1100, 700

    def __init__(self, page: "FigurePage", key: tuple[str, str], opts, current=None):
        from ..cellfigure import condition_fields

        super().__init__(page)
        self.page, self.key, self.opts = page, key, opts
        self.cd = page.cd
        self.fields = condition_fields(self.cd, *key)
        self.idx = 0
        if current is not None:
            self.idx = next((i for i, f in enumerate(self.fields) if (f.well, f.field) == (current.well, current.field)), 0)
        self.sel = current.nucleus_id if current is not None else 0
        self.factor, self.scale, self.small_mask, self.base = 4, 1.0, None, None
        self.title(f"Choose a cell · {key[0]} · {key[1]}")
        self.geometry("1180x900")
        self.configure(fg_color=T.BG)
        self.transient(page.winfo_toplevel())

        top = ctk.CTkFrame(self, fg_color="transparent")
        top.pack(fill="x", padx=16, pady=(12, 6))
        ghost_button(top, "◀", lambda: self.go(-1), width=40).pack(side="left")
        self.field_menu = _menu(top, [self._flabel(f) for f in self.fields] or ["(no fields)"], 340,
                                command=self._field_chosen)
        self.field_menu.pack(side="left", padx=6)
        ghost_button(top, "▶", lambda: self.go(1), width=40).pack(side="left")
        signal = [c for c in self.cd.ds.channels if c != self.cd.nuclear]
        self.views = {"Merge + nuclear stain": (*signal, self.cd.nuclear) if self.cd.nuclear else tuple(signal),
                      "Merge": tuple(signal), **{c: (c,) for c in self.cd.ds.channels}}
        self.view = _menu(top, list(self.views), 200, command=lambda v: self.load_field())
        self.view.pack(side="left", padx=(16, 0))
        ctk.CTkLabel(top, text="white = analysed · red = QC-flagged · yellow = chosen · box = crop · #1 = suggestions",
                     font=T.font(11), text_color=T.MUTED).pack(side="right")

        import tkinter as tk

        self.canvas = tk.Canvas(self, width=self.MAX_W, height=self.MAX_H, highlightthickness=0, bd=0, bg="#000000",
                                cursor="crosshair")
        self.canvas.pack(padx=16)
        self.canvas.bind("<Button-1>", self.click)

        bottom = ctk.CTkFrame(self, fg_color="transparent")
        bottom.pack(fill="x", padx=16, pady=10)
        self.info = ctk.CTkLabel(bottom, text="Click a nucleus.", font=T.font(12), text_color=T.TEXT, anchor="w",
                                 justify="left", wraplength=800)
        self.info.pack(side="left", fill="x", expand=True)
        ghost_button(bottom, "Cancel", self.destroy, width=90).pack(side="right")
        self.use_btn = primary_button(bottom, "Use this cell", self.use, width=130, state="disabled")
        self.use_btn.pack(side="right", padx=8)
        if self.fields:
            self.field_menu.set(self._flabel(self.fields[self.idx]))
            self.load_field()
        self.after(100, self.lift)

    def _flabel(self, f) -> str:
        cells = self.cd.cells
        n = int(((cells["well"] == f.well) & (cells["field"] == f.field)).sum())
        return f"{f.well} · {f.field}  ({n} analysed nuclei)"

    def _field_chosen(self, label: str) -> None:
        self.idx = [self._flabel(f) for f in self.fields].index(label)
        self.sel = 0
        self.load_field()

    def go(self, step: int) -> None:
        if not self.fields:
            return
        self.idx = (self.idx + step) % len(self.fields)
        self.sel = 0
        self.field_menu.set(self._flabel(self.fields[self.idx]))
        self.load_field()

    @property
    def field(self):
        return self.fields[self.idx]

    def load_field(self) -> None:
        from ..cellfigure import field_overview

        f, view = self.field, self.views[self.view.get()]
        self.info.configure(text=f"Loading {f.well} {f.field}…", text_color=T.MUTED)
        colors = dict(self.opts.colors)
        if self.cd.nuclear:
            colors.setdefault(self.cd.nuclear, "Blue")
            if colors[self.cd.nuclear] == "Gray" and len(view) > 1:
                colors[self.cd.nuclear] = "Blue"  # gray nuclei wash out the merge in the overview

        def work():
            with self.page.lock:
                return field_overview(self.cd, f.well, f.field, view, colors, self.factor, self.page.cache)

        def done(res):
            if not self.winfo_exists() or self.field is not f:
                return
            self.base, self.small_mask = res
            h, w = self.small_mask.shape
            self.scale = min(self.MAX_W / w, self.MAX_H / h, 2.0)
            self.canvas.configure(width=int(w * self.scale), height=int(h * self.scale))
            self.redraw()

        self.page.app.run_bg(work, done, lambda exc: self.info.configure(text=f"Could not load: {exc}",
                                                                        text_color=T.DANGER))

    def redraw(self) -> None:
        from PIL import Image, ImageTk
        from skimage.segmentation import find_boundaries

        if self.base is None:
            return
        f = self.field
        cells = self.cd.cells
        ok = set(cells.loc[(cells["well"] == f.well) & (cells["field"] == f.field), "nucleus_id"].astype(int))
        m = self.small_mask
        edges = find_boundaries(m, mode="inner")
        rgb = self.base.copy()
        labels = m[edges]
        good = np.isin(labels, list(ok)) if ok else np.zeros(labels.shape, bool)
        col = np.where(good[:, None], np.array([235, 235, 235]), np.array([230, 60, 60])).astype(np.uint8)
        rgb[edges] = col
        if self.sel:
            from scipy import ndimage as ndi

            own = m == self.sel
            rgb[own & ~ndi.binary_erosion(own, iterations=2)] = (255, 212, 0)
        h, w = m.shape
        img = Image.fromarray(rgb).resize((int(w * self.scale), int(h * self.scale)), Image.NEAREST
                                          if self.scale >= 1 else Image.BILINEAR)
        self._photo = ImageTk.PhotoImage(img)
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, image=self._photo, anchor="nw")
        # ranked suggestions in this field, labelled #1, #2 ... (#1 = automatic pick)
        cands = self.page.candidates.get(self.key)
        if cands is not None and len(cands):
            px = self.cd.ds.px_um
            for i, (_, r) in enumerate(cands.iterrows()):
                if (r["well"], r["field"]) == (f.well, f.field) and str(r.get("relax", "")) != "nan":
                    x = r["centroid_x_um"] / px / self.factor * self.scale
                    y = r["centroid_y_um"] / px / self.factor * self.scale
                    for dx, dy, color in ((1, 1, "#000000"), (0, 0, "#7CF2FF")):
                        self.canvas.create_text(x + dx, y + dy - 26 * self.scale, text=f"#{i + 1}", fill=color,
                                                font=("Segoe UI", 10, "bold"))
        if self.sel:
            self._draw_crop()
        self._update_info()

    def _draw_crop(self) -> None:
        own = np.nonzero(self.small_mask == self.sel)
        if not len(own[0]):
            return
        px = self.cd.ds.px_um
        cy, cx = float(own[0].mean()), float(own[1].mean())
        hh = self.opts.crop_h_um / px / self.factor / 2
        hw = self.opts.crop_w_um / px / self.factor / 2
        s = self.scale
        self.canvas.create_rectangle((cx - hw) * s, (cy - hh) * s, (cx + hw) * s, (cy + hh) * s, outline="#FFD400",
                                     width=2, dash=(6, 4))

    def click(self, event) -> None:
        from ..cellfigure import nucleus_at

        if self.small_mask is None:
            return
        y = int(event.y / self.scale * self.factor)
        x = int(event.x / self.scale * self.factor)
        f = self.field
        nid = nucleus_at(self.cd, f.well, f.field, y, x, radius=int(2 / self.cd.ds.px_um), cache=self.page.cache)
        self.sel = nid
        self.redraw()

    def _update_info(self) -> None:
        from ..cellfigure import cell_summary, describe

        f = self.field
        if not self.sel:
            self.info.configure(text="Click a nucleus.", text_color=T.MUTED)
            self.use_btn.configure(state="disabled")
            return
        r = cell_summary(self.cd, *self.key, f.well, f.field, self.sel)
        if r is None:
            self.info.configure(text=f"{f.well} {f.field} #{self.sel}: flagged by QC (condensed, dim or saturated) or "
                                     "not analysed - choose another nucleus.", text_color=T.WARNING)
            self.use_btn.configure(state="disabled")
            return
        warn = []
        if r.get("dapi_ratio", 1.0) > 1.25 or r.get("dapi_ratio", 1.0) < 0.75:
            warn.append("DAPI suggests a condensed / dim nucleus")
        px = self.cd.ds.px_um
        cy, cx = r["centroid_y_um"] / px, r["centroid_x_um"] / px
        H, W = self.small_mask.shape[0] * self.factor, self.small_mask.shape[1] * self.factor
        if (cy - self.opts.crop_h_um / px / 2 < 0 or cx - self.opts.crop_w_um / px / 2 < 0
                or cy + self.opts.crop_h_um / px / 2 > H or cx + self.opts.crop_w_um / px / 2 > W):
            warn.append("the crop runs off the image (padded black)")
        self.info.configure(text="cell (condition median):  " + describe(r, self.cd) +
                                 (("\n⚠ " + "; ".join(warn)) if warn else ""),
                            text_color=T.WARNING if warn else T.TEXT)
        self._row = r
        self.use_btn.configure(state="normal")

    def use(self) -> None:
        self.page.set_custom(self.key, self._row)
        self.destroy()


class FigurePage(ctk.CTkFrame):
    title = "Figures"
    subtitle = "Publication single-cell figures - one representative cell per cell type and treatment"

    def __init__(self, master, app):
        super().__init__(master, fg_color="transparent")
        from ..cellfigure import PAGE_WIDTHS_IN

        self.app = app
        self.cd = None
        self.out: "Path | None" = None
        self.candidates: dict = {}  # (ct, t) -> DataFrame
        self.rows_ui: dict = {}  # (ct, t) -> (menu, thumb label, info label)
        self.cache = None
        self.lock = threading.Lock()
        self.layout_card = None
        self.body = ctk.CTkScrollableFrame(self, fg_color="transparent", scrollbar_button_color=T.BORDER,
                                           scrollbar_button_hover_color=T.MUTED)
        self.body.pack(fill="both", expand=True)

        # input
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Input", "A Zeizify output folder that has been analysed (Analyze page)").pack(
            fill="x", padx=18, pady=(14, 6))
        r = ctk.CTkFrame(c, fg_color="transparent")
        r.pack(fill="x", padx=18, pady=(0, 6))
        self.folder_btn = ctk.CTkButton(r, text="▢  Choose output folder…", anchor="w", height=32, corner_radius=16,
                                        fg_color=T.SURFACE_ALT, hover_color=T.BORDER, text_color=T.MUTED,
                                        font=T.font(12), command=self.pick_folder)
        self.folder_btn.pack(side="left", fill="x", expand=True)
        self.info = ctk.CTkLabel(c, text="", font=T.font(12), text_color=T.MUTED, anchor="w", justify="left",
                                 wraplength=900)
        self.info.pack(fill="x", padx=18, pady=(0, 14))

        # cells
        c = self.cells_card = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Cells", "Conditions to show, the crop size, and one cell per condition. Pick "
                                  "representative cells, then change any of them from its ranked list.").pack(
            fill="x", padx=18, pady=(14, 6))
        self.ct_row = _row(c, "Cell types")
        self.tr_row = _row(c, "Treatments")
        self.ct_checks: dict[str, ctk.CTkCheckBox] = {}
        self.tr_checks: dict[str, ctk.CTkCheckBox] = {}
        r = _row(c, "Cell-type order / labels", "order: comma-separated · labels: WT=Wild type; KO1=PML KO")
        self.ct_order, self.ct_labels = _entry(r, "", 220), _entry(r, "", 300)
        self.ct_order.pack(side="left")
        self.ct_labels.pack(side="left", padx=8)
        r = _row(c, "Treatment order / labels", "e.g. ML=ML792; MG=MG132; ML+MG=ML792 + MG132")
        self.tr_order, self.tr_labels = _entry(r, "", 220), _entry(r, "", 300)
        self.tr_order.pack(side="left")
        self.tr_labels.pack(side="left", padx=8)
        r = _row(c, "Crop size (µm)", "width × height, centred on the nucleus")
        self.crop_w, self.crop_h = _entry(r, 30, 60), _entry(r, 30, 60)
        self.crop_w.pack(side="left")
        ctk.CTkLabel(r, text="×", text_color=T.MUTED).pack(side="left", padx=6)
        self.crop_h.pack(side="left")
        r = ctk.CTkFrame(c, fg_color="transparent")
        r.pack(fill="x", padx=18, pady=(8, 4))
        self.pick_btn = primary_button(r, "Pick representative cells", self.pick_cells, width=200)
        self.pick_btn.pack(side="left")
        self.pick_status = ctk.CTkLabel(r, text="", font=T.font(11), text_color=T.MUTED)
        self.pick_status.pack(side="left", padx=10)
        self.table = ctk.CTkFrame(c, fg_color="transparent")
        self.table.pack(fill="x", padx=18, pady=(4, 14))

        # panels + colors (the Contrast page's layout card) and their order
        self.layout_slot = ctk.CTkFrame(self.body, fg_color="transparent")
        self.layout_slot.pack(fill="x")
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Layout", "Panel order, orientation and page size").pack(fill="x", padx=18, pady=(14, 6))
        r = _row(c, "Panel order", "edit to reorder: DAPI | AF488 | AF488+AF594")
        self.order = _entry(r, "", 420)
        self.order.pack(side="left")
        r = _row(c, "Arrangement")
        self.rows = _seg(r, ROWS)
        self.rows.pack(side="left")
        self.blocks = _seg(r, BLOCKS)
        self.blocks.pack(side="left", padx=10)
        r = _row(c, "Figure width")
        self.width_menu = _menu(r, [*PAGE_WIDTHS_IN, CUSTOM_WIDTH], 260)
        self.width_menu.set("Journal - 2 columns (183 mm)")
        self.width_menu.pack(side="left")
        self.width_mm = _entry(r, "", 70)
        self.width_mm.pack(side="left", padx=8)
        ctk.CTkLabel(r, text="mm (custom)", font=T.font(11), text_color=T.MUTED).pack(side="left")
        r = _row(c, "Font", "Arial is accepted by most journals")
        self.font = _menu(r, FONTS, 160)
        self.font.pack(side="left")
        self.font_size = _entry(r, 8, 50)
        self.font_size.pack(side="left", padx=8)
        ctk.CTkLabel(r, text="pt", font=T.font(11), text_color=T.MUTED).pack(side="left")
        r = _row(c, "Gaps (mm)", "between panels · between cell-type blocks")
        self.gap, self.block_gap = _entry(r, 0.5, 50), _entry(r, 3, 50)
        self.gap.pack(side="left")
        self.block_gap.pack(side="left", padx=8)
        ctk.CTkFrame(c, fg_color="transparent", height=8).pack()

        # display
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Display", "Linear, background-subtracted. Use the same range in every condition so "
                                    "intensities can be compared between panels.").pack(fill="x", padx=18, pady=(14, 6))
        r = _row(c, "Brightness range")
        self.contrast = _seg(r, CONTRAST)
        self.contrast.pack(side="left")
        r = _row(c, "Saturated pixels (%)", "share of nuclear pixels allowed to saturate (sets the maximum)")
        self.sat = _entry(r, 0.3, 60)
        self.sat.pack(side="left")
        self.max_row = _row(c, "Maximum per channel", "blank = automatic")
        self.max_entries: dict[str, ctk.CTkEntry] = {}
        r = _row(c, "Outline measured nucleus")
        self.outline = _menu(r, OUTLINE, 200)
        self.outline.pack(side="left")
        r = _row(c, "Scale bar")
        self.sb_len = _entry(r, 5, 50)
        self.sb_len.pack(side="left")
        ctk.CTkLabel(r, text="µm", font=T.font(11), text_color=T.MUTED).pack(side="left", padx=(4, 10))
        self.sb_where = _menu(r, SCALEBAR, 180)
        self.sb_where.pack(side="left")
        self.sb_label = _check(r, "Show length on the bar", True)
        self.sb_label.pack(side="left", padx=12)
        ctk.CTkFrame(c, fg_color="transparent", height=8).pack()

        # export
        c = Card(self.body)
        c.pack(fill="x", pady=(0, 12))
        section_title(c, "Export", "Writes <output folder>/figures/<name>/: the figure, cells.csv, legend.txt and "
                                   "figure_settings.json").pack(fill="x", padx=18, pady=(14, 6))
        r = _row(c, "Name")
        self.name = _entry(r, "single_cells", 260)
        self.name.pack(side="left")
        r = _row(c, "Formats", "PDF / SVG keep text editable")
        self.fmt_checks = {f: _check(r, f.upper(), f in ("png", "pdf")) for f in ("png", "tif", "pdf", "svg")}
        for cb in self.fmt_checks.values():
            cb.pack(side="left", padx=(0, 12))
        r = _row(c, "Resolution")
        self.dpi = _menu(r, ["300 dpi", "600 dpi"], 120)
        self.dpi.pack(side="left")
        r = ctk.CTkFrame(c, fg_color="transparent")
        r.pack(fill="x", padx=18, pady=(8, 6))
        self.preview_btn = ghost_button(r, "Preview", self.preview, width=100)
        self.preview_btn.pack(side="left")
        self.export_btn = primary_button(r, "Export figure", self.export, width=140)
        self.export_btn.pack(side="left", padx=8)
        self.open_btn = ghost_button(r, "Open folder", self.open_folder, width=110)
        self.open_btn.pack(side="left")
        self.status = ctk.CTkLabel(c, text="", font=T.font(12), text_color=T.MUTED, anchor="w", justify="left",
                                   wraplength=900)
        self.status.pack(fill="x", padx=18)
        self.preview_img = ctk.CTkLabel(c, text="")
        self.preview_img.pack(padx=18, pady=(6, 14))
        self._last_folder = None

    # ---------- input ----------

    def pick_folder(self) -> None:
        d = filedialog.askdirectory(title="Zeizify output folder (analysed)")
        if d:
            self.load(Path(d))

    def on_show(self) -> None:
        if self.cd is not None:
            return
        for page in ("analyze", "run"):
            out = getattr(self.app.pages.get(page), "out", None)
            if out and (Path(out) / "analysis" / "nuclei.csv").exists():
                self.load(Path(out))
                return

    def load(self, out: Path) -> None:
        from ..analysis.stats import is_control, is_wildtype, level_order
        from ..cellfigure import _FieldCache, default_figure_panels, load_cells
        from .contrast_page import FigureLayoutCard

        try:
            cd = load_cells(out)
        except Exception as exc:  # noqa: BLE001
            self.info.configure(text=f"Could not load: {exc}", text_color=T.DANGER)
            return
        self.cd, self.out, self.candidates = cd, out, {}
        self.cache = _FieldCache(cd, size=6)
        self.folder_btn.configure(text=f"▢  {out}", text_color=T.TEXT)
        n = len(cd.cells)
        self.info.configure(text=f"{n} analysed nuclei (QC-passed) · {len(cd.ds.fields)} fields · channels: "
                                 f"{', '.join(cd.ds.channels)} · {cd.ds.px_um:.4f} µm/px", text_color=T.MUTED)
        cts = level_order(cd.cells, "cell_type", first=is_wildtype)
        trs = level_order(cd.cells, "treatment", first=is_control)
        for row, checks, values in ((self.ct_row, self.ct_checks, cts), (self.tr_row, self.tr_checks, trs)):
            for w in row.winfo_children()[1:]:
                w.destroy()
            checks.clear()
            for v in values:
                cb = _check(row, v or "(none)", True)
                cb.pack(side="left", padx=(0, 12))
                checks[v] = cb
        for e, v in ((self.ct_order, ", ".join(cts)), (self.tr_order, ", ".join(trs))):
            e.delete(0, "end")
            e.insert(0, v)
        for w in self.table.winfo_children():
            w.destroy()
        self.rows_ui = {}
        self.pick_status.configure(text="")
        # panels / colors
        for w in self.layout_slot.winfo_children():
            w.destroy()
        self.layout_card = FigureLayoutCard(
            self.layout_slot, cd.ds.channels, self._panels_changed, title="Panels and colors",
            subtitle="Tick the panels to show. Merges add the channels' colors (green + magenta = white where they "
                     "overlap). The order is set under Layout.", defaults=default_figure_panels)
        self.layout_card.pack(fill="x", pady=(0, 12))
        self.layout_card._preset(["Green", "Magenta", "Cyan", "Yellow"])
        if cd.nuclear:
            self.layout_card.color_menus[cd.nuclear].set("Gray")
        self._set_order(default_figure_panels(cd.ds.channels))
        for w in self.max_row.winfo_children()[1:]:
            w.destroy()
        self.max_entries = {}
        for ch in cd.ds.channels:
            box = ctk.CTkFrame(self.max_row, fg_color=T.SURFACE_ALT, corner_radius=8)
            box.pack(side="left", padx=(0, 8))
            ctk.CTkLabel(box, text=ch, font=T.font(11), text_color=T.MUTED).pack(side="left", padx=(8, 4))
            e = ctk.CTkEntry(box, width=70, height=26, corner_radius=6, border_width=0, fg_color=T.SURFACE,
                             font=T.font(11))
            e.pack(side="left", padx=(0, 6), pady=4)
            self.max_entries[ch] = e

    def _set_order(self, panels) -> None:
        self.order.delete(0, "end")
        self.order.insert(0, " | ".join("+".join(p) for p in panels))

    def _panels_changed(self) -> None:
        """Keep the order of panels that stay ticked; append new ones at the end."""
        from ..cellfigure import parse_order

        ticked = self.layout_card.panels()
        cur = [p for p in parse_order(self.order.get()) if p in ticked]
        self._set_order(cur + [p for p in ticked if p not in cur])

    # ---------- options ----------

    def _ordered(self, entry, checks) -> list[str]:
        from ..analysis.stats import _ordered

        typed = [t.strip() for t in entry.get().split(",") if t.strip()]
        on = {k for k, cb in checks.items() if cb.get()}
        return [x for x in _ordered(on, typed)]

    def options(self):
        from ..cellfigure import PAGE_WIDTHS_IN, FigureOptions, parse_labels, parse_order

        try:
            panels = parse_order(self.order.get())
            bad = [c for p in panels for c in p if c not in self.cd.ds.channels]
            if bad:
                raise ValueError(f"unknown channel(s) in the panel order: {', '.join(bad)}")
            if not panels:
                raise ValueError("no panels - tick at least one under Panels and colors")
            w = self.width_menu.get()
            width_in = float(self.width_mm.get()) / 25.4 if w == CUSTOM_WIDTH else PAGE_WIDTHS_IN[w]
            maxes = {}
            for ch, e in self.max_entries.items():
                if e.get().strip():
                    maxes[ch] = float(e.get())
            return FigureOptions(
                panels=panels, colors=self.layout_card.colors(),
                cell_types=self._ordered(self.ct_order, self.ct_checks),
                treatments=self._ordered(self.tr_order, self.tr_checks),
                cell_type_labels=parse_labels(self.ct_labels.get()), treatment_labels=parse_labels(self.tr_labels.get()),
                crop_w_um=float(self.crop_w.get()), crop_h_um=float(self.crop_h.get()),
                rows=ROWS[self.rows.get()], blocks=BLOCKS[self.blocks.get()], width_in=width_in,
                dpi=int(self.dpi.get().split()[0]), font=self.font.get(), font_size=float(self.font_size.get()),
                gap_mm=float(self.gap.get()), block_gap_mm=float(self.block_gap.get()),
                contrast=CONTRAST[self.contrast.get()], saturated_pct=float(self.sat.get()), max_override=maxes,
                outline=OUTLINE[self.outline.get()], scale_bar_um=float(self.sb_len.get()),
                scale_bar=SCALEBAR[self.sb_where.get()], scale_bar_label=bool(self.sb_label.get()),
                formats=[f for f, cb in self.fmt_checks.items() if cb.get()] or ["png"])
        except ValueError as exc:
            self.app.toast(f"Check the settings: {exc}", kind="warning")
            return None

    # ---------- cells ----------

    def pick_cells(self) -> None:
        from ..cellfigure import auto_choices

        if self.cd is None:
            self.app.toast("Choose an analysed output folder first.", kind="warning")
            return
        opts = self.options()
        if opts is None:
            return
        conds = [(c, t) for c in opts.cell_types for t in opts.treatments if (c, t) in set(self.cd.conditions())]
        if not conds:
            self.app.toast("No condition selected.", kind="warning")
            return
        self.pick_btn.configure(state="disabled", text="Ranking…")

        def progress(i, n, item):
            self.app.post(lambda: self.pick_status.configure(text=f"Ranking cells {i} / {n} · {item}"))

        def work():
            with self.lock:
                return auto_choices(self.cd, conds, opts.crop_w_um, opts.crop_h_um, n=12, progress=progress)

        def done(res):
            self.candidates = res
            self.pick_btn.configure(state="normal", text="Pick representative cells")
            empty = [f"{c} · {t}" for (c, t), v in res.items() if not len(v)]
            self.pick_status.configure(text=f"{len(res) - len(empty)} conditions ready" +
                                       (f" · no suitable cell: {', '.join(empty)}" if empty else ""),
                                       text_color=T.WARNING if empty else T.MUTED)
            self._build_table(opts)

        def fail(exc):
            self.pick_btn.configure(state="normal", text="Pick representative cells")
            self.pick_status.configure(text=f"Failed: {exc}", text_color=T.DANGER)

        self.app.run_bg(work, done, fail)

    def _build_table(self, opts) -> None:
        from ..cellfigure import describe

        for w in self.table.winfo_children():
            w.destroy()
        self.rows_ui = {}
        for (ct, t), cands in self.candidates.items():
            if not len(cands):
                continue
            row = ctk.CTkFrame(self.table, fg_color=T.SURFACE_ALT, corner_radius=10)
            row.pack(fill="x", pady=3)
            thumb = ctk.CTkLabel(row, text="", width=72, height=72)
            thumb.pack(side="left", padx=8, pady=6)
            col = ctk.CTkFrame(row, fg_color="transparent")
            col.pack(side="left", fill="x", expand=True, padx=(0, 8))
            ctk.CTkLabel(col, text=f"{ct} · {t}", font=T.font(12, "bold"), text_color=T.TEXT, anchor="w").pack(fill="x")
            labels = [f"#{i + 1}{' (auto)' if i == 0 else ''}  {describe(r, self.cd)}" for i, (_, r) in
                      enumerate(cands.iterrows())]
            line = ctk.CTkFrame(col, fg_color="transparent")
            line.pack(fill="x", pady=(2, 0))
            menu = ctk.CTkOptionMenu(line, values=labels, height=28, corner_radius=8, font=T.font(11),
                                     fg_color=T.SURFACE, button_color=T.SURFACE, button_hover_color=T.BORDER,
                                     text_color=T.TEXT, dropdown_font=T.font(11), dynamic_resizing=False, width=640,
                                     command=lambda v, k=(ct, t): self._chosen(k))
            menu.set(labels[0])
            menu.pack(side="left")
            ghost_button(line, "Choose on image…", lambda k=(ct, t): self.open_picker(k), width=140,
                         height=28).pack(side="left", padx=8)
            note = ctk.CTkLabel(col, text="cell (condition median)", font=T.font(10), text_color=T.MUTED, anchor="w")
            note.pack(fill="x")
            self.rows_ui[(ct, t)] = (menu, thumb, labels)
            self._chosen((ct, t))

    def _choice_index(self, k) -> int:
        menu, _, labels = self.rows_ui[k]
        return labels.index(menu.get()) if menu.get() in labels else 0

    def _chosen(self, k) -> None:
        from ..cellfigure import choice_from_row, thumbnail

        opts = self.options()
        if opts is None:
            return
        i = self._choice_index(k)
        choice = choice_from_row(self.candidates[k].iloc[i], "auto" if i == 0 else "manual")
        _, thumb, _ = self.rows_ui[k]

        def work():
            with self.lock:
                return thumbnail(self.cd, choice, opts, 72, self.cache)

        def done(img):
            if k in self.rows_ui and self.rows_ui[k][1].winfo_exists():
                self.rows_ui[k][1].configure(image=ctk.CTkImage(img, img, size=img.size))

        self.app.run_bg(work, done, lambda exc: None)

    def open_picker(self, k) -> None:
        opts = self.options()
        if opts is None:
            return
        if getattr(self, "_picker", None) is not None and self._picker.winfo_exists():
            self._picker.destroy()
        self._picker = FieldPicker(self, k, opts, self.choices().get(k))

    def set_custom(self, k, row) -> None:
        """Use a nucleus chosen on the field image: add it to the condition's list and select it."""
        import pandas as pd

        from ..cellfigure import describe

        cands = self.candidates[k]
        ident = (row["well"], row["field"], int(row["nucleus_id"]))
        hit = [i for i, (_, r) in enumerate(cands.iterrows()) if (r["well"], r["field"], int(r["nucleus_id"])) == ident]
        if hit:
            i = hit[0]
        else:
            self.candidates[k] = pd.concat([cands, row.to_frame().T], ignore_index=True)
            i = len(self.candidates[k]) - 1
        menu, thumb, labels = self.rows_ui[k]
        label = f"#{i + 1} (chosen on image)  {describe(row, self.cd)}"
        labels = labels[:i] + [label] + labels[i + 1:] if i < len(labels) else labels + [label]
        menu.configure(values=labels)
        menu.set(label)
        self.rows_ui[k] = (menu, thumb, labels)
        self._chosen(k)
        self.app.toast(f"{k[0]} · {k[1]}: using {ident[0]} {ident[1]} #{ident[2]}", kind="success")

    def choices(self) -> dict:
        from ..cellfigure import choice_from_row

        out = {}
        for k in self.rows_ui:
            i = self._choice_index(k)
            out[k] = choice_from_row(self.candidates[k].iloc[i], "auto" if i == 0 else "manual")
        return out

    # ---------- preview / export ----------

    def _ready(self):
        if self.cd is None:
            self.app.toast("Choose an analysed output folder first.", kind="warning")
            return None
        if not self.rows_ui:
            self.app.toast("Press 'Pick representative cells' first.", kind="warning")
            return None
        return self.options()

    def preview(self) -> None:
        from ..cellfigure import render

        opts = self._ready()
        if opts is None:
            return
        choices = self.choices()
        self.preview_btn.configure(state="disabled", text="Drawing…")
        self.status.configure(text="Drawing the preview…", text_color=T.MUTED)

        def work():
            import matplotlib.pyplot as plt

            with self.lock:
                fig, _, _ = render(self.cd, choices, opts)
            buf = io.BytesIO()
            fig.savefig(buf, dpi=110, facecolor="white")
            w_in, h_in = fig.get_size_inches()
            plt.close(fig)
            buf.seek(0)
            return Image.open(buf), (w_in, h_in)

        def done(res):
            img, (w_in, h_in) = res
            self.preview_btn.configure(state="normal", text="Preview")
            w = min(1000, img.width)
            self.preview_img.configure(image=ctk.CTkImage(img, img, size=(w, int(img.height * w / img.width))))
            self.status.configure(text=f"Preview · {w_in * 25.4:.0f} × {h_in * 25.4:.0f} mm "
                                       f"({w_in:.2f} × {h_in:.2f} in) at {opts.dpi} dpi on export", text_color=T.MUTED)

        def fail(exc):
            self.preview_btn.configure(state="normal", text="Preview")
            self.status.configure(text=f"Preview failed: {exc}", text_color=T.DANGER)

        self.app.run_bg(work, done, fail)

    def export(self) -> None:
        from ..cellfigure import export

        opts = self._ready()
        if opts is None:
            return
        choices = self.choices()
        name = self.name.get().strip() or "single_cells"
        self.export_btn.configure(state="disabled", text="Exporting…")

        def work():
            with self.lock:
                return export(self.cd, choices, opts, name, self.candidates)

        def done(res):
            self.export_btn.configure(state="normal", text="Export figure")
            self._last_folder = res["folder"]
            self.status.configure(text=f"Saved to {res['folder']}:  " + ", ".join(p.name for p in res["files"]) +
                                       ", cells.csv, legend.txt, figure_settings.json", text_color=T.SUCCESS)
            self.app.toast("Figure exported.", kind="success")

        def fail(exc):
            self.export_btn.configure(state="normal", text="Export figure")
            self.status.configure(text=f"Export failed: {exc}", text_color=T.DANGER)

        self.app.run_bg(work, done, fail)

    def open_folder(self) -> None:
        from .run_page import open_path

        target = self._last_folder or (self.out / "figures" if self.out else None)
        if target and Path(target).exists():
            open_path(Path(target))
        else:
            self.app.toast("Nothing exported yet.", kind="info")

    def on_theme_change(self) -> None:
        pass
