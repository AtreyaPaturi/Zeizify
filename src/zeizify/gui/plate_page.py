"""
Plate map editor: paint wells with three independent categories - A: cell type,
B: treatment, C: antibody - and export the CSV the extraction step reads.

Pick a label on the right, then click or drag across wells to paint (clicking a well
that already has that label clears it). Click a row letter, column number or the
corner to toggle a whole row/column/plate. Each well shows its three categories as
three wedges.
"""

from __future__ import annotations

import tkinter as tk
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk

from ..plate import PLATE_FORMATS
from . import theme as T
from .state import CATEGORIES, CATEGORY_BADGES, CATEGORY_TITLES, PlateState
from .widgets import Card, ghost_button, primary_button, section_title

WEDGE_START = {"cellType": 90, "treatment": 210, "antibody": 330}


class PlateCanvas(tk.Canvas):
    def __init__(self, master, state: PlateState, on_change, on_hover):
        super().__init__(master, highlightthickness=0, bd=0)
        self.state = state
        self.on_change, self.on_hover_cb = on_change, on_hover
        self.brush: "tuple[str, str] | None" = None  # (category, label id | "erase")
        self.drag_mode: "str | None" = None
        self.wedges: dict[str, dict[str, int]] = {}
        self.bind("<ButtonPress-1>", self._press)
        self.bind("<B1-Motion>", self._drag)
        self.bind("<ButtonRelease-1>", lambda e: setattr(self, "drag_mode", None))
        self.bind("<Motion>", self._hover)
        self.bind("<Leave>", lambda e: self.on_hover_cb(None))
        self.redraw()

    # geometry
    @property
    def cell(self) -> int:
        return max(20, min(46, 600 // self.state.plate.n_cols))

    @property
    def header(self) -> int:
        return 28

    def redraw(self) -> None:
        self.delete("all")
        p, cell, hd = self.state.plate, self.cell, self.header
        self.configure(width=hd + p.n_cols * cell + 4, height=hd + p.n_rows * cell + 4, bg=T.pick(T.SURFACE))
        small = (T.FAMILY, 8 if cell < 30 else 9)
        self.create_text(hd / 2, hd / 2, text="all", fill=T.pick(T.MUTED), font=small)
        for c in range(1, p.n_cols + 1):
            self.create_text(hd + (c - 0.5) * cell, hd / 2, text=str(c), fill=T.pick(T.MUTED), font=small)
        self.wedges.clear()
        pad = max(2, cell // 10)
        for ri, r in enumerate(p.row_letters):
            y0 = hd + ri * cell
            self.create_text(hd / 2, y0 + cell / 2, text=r, fill=T.pick(T.MUTED), font=small)
            for c in range(1, p.n_cols + 1):
                x0 = hd + (c - 1) * cell
                box = (x0 + pad, y0 + pad, x0 + cell - pad, y0 + cell - pad)
                well = f"{r}{c}"
                self.wedges[well] = {
                    cat: self.create_arc(*box, start=WEDGE_START[cat], extent=120, style="pieslice",
                                         outline=T.pick(T.SURFACE), width=1)
                    for cat in CATEGORIES
                }
                self.create_oval(*box, outline=T.pick(T.WELL_LINE), width=1)
        self.paint_all()

    def paint(self, well: str) -> None:
        for cat, item in self.wedges[well].items():
            self.itemconfigure(item, fill=self.state.color_at(cat, well) or T.pick(T.WELL_EMPTY))

    def paint_all(self) -> None:
        for w in self.wedges:
            self.paint(w)

    # hit testing
    def _well_at(self, x, y) -> "str | None":
        p, cell, hd = self.state.plate, self.cell, self.header
        if x < hd or y < hd:
            return None
        c, r = int((x - hd) // cell), int((y - hd) // cell)
        return f"{p.row_letters[r]}{c + 1}" if 0 <= r < p.n_rows and 0 <= c < p.n_cols else None

    def _header_wells(self, x, y) -> "list[str] | None":
        p, cell, hd = self.state.plate, self.cell, self.header
        if x < hd and y < hd:
            return p.wells()
        if y < hd:
            c = int((x - hd) // cell) + 1
            return p.col_wells(c) if 1 <= c <= p.n_cols else None
        if x < hd:
            r = int((y - hd) // cell)
            return p.row_wells(p.row_letters[r]) if 0 <= r < p.n_rows else None
        return None

    # painting
    def _apply(self, well: str, mode: str) -> None:
        cat, label = self.brush
        by_well = self.state.cats[cat]["byWell"]
        if mode == "erase":
            by_well.pop(well, None)
        else:
            by_well[well] = label
        self.paint(well)

    def _press(self, e) -> None:
        if not self.brush:
            self.on_hover_cb("Pick a label on the right first, then paint wells.")
            return
        cat, label = self.brush
        by_well = self.state.cats[cat]["byWell"]
        line = self._header_wells(e.x, e.y)
        if line:
            all_set = label != "erase" and all(by_well.get(w) == label for w in line)
            for w in line:
                self._apply(w, "erase" if (label == "erase" or all_set) else "paint")
            self.on_change(cat)
            return
        well = self._well_at(e.x, e.y)
        if well:
            self.drag_mode = "erase" if label == "erase" or by_well.get(well) == label else "paint"
            self._apply(well, self.drag_mode)
            self.on_change(cat)

    def _drag(self, e) -> None:
        well = self._well_at(e.x, e.y)
        if self.drag_mode and well:
            self._apply(well, self.drag_mode)
            self.on_change(self.brush[0])

    def _hover(self, e) -> None:
        well = self._well_at(e.x, e.y)
        if not well:
            self.on_hover_cb(None)
            return
        parts = [f"{CATEGORY_TITLES[c]}: {self.state.name_at(c, well) or '—'}" for c in CATEGORIES]
        self.on_hover_cb(f"{well}   ·   " + "   ·   ".join(parts))


class CategoryPanel(Card):
    def __init__(self, master, page: "PlatePage", cat: str):
        super().__init__(master)
        self.page, self.cat = page, cat
        head = ctk.CTkFrame(self, fg_color="transparent")
        head.pack(fill="x", padx=16, pady=(14, 4))
        ctk.CTkLabel(head, text=CATEGORY_BADGES[cat], width=24, height=24, corner_radius=6,
                     fg_color=T.ACCENT_SOFT, text_color=T.ACCENT, font=T.font(12, "bold")).pack(side="left")
        ctk.CTkLabel(head, text=CATEGORY_TITLES[cat], font=T.font(14, "bold"), text_color=T.TEXT).pack(side="left", padx=8)
        self.sub = ctk.CTkLabel(head, text="", font=T.font(11), text_color=T.MUTED)
        self.sub.pack(side="right")

        self.chips = ctk.CTkFrame(self, fg_color="transparent")
        self.chips.pack(fill="x", padx=12)

        add = ctk.CTkFrame(self, fg_color="transparent")
        add.pack(fill="x", padx=16, pady=(6, 14))
        self.entry = ctk.CTkEntry(add, placeholder_text=f"New {CATEGORY_TITLES[cat].lower()} label…", height=32,
                                  corner_radius=8, border_color=T.BORDER, fg_color=T.SURFACE_ALT, font=T.font(12))
        self.entry.pack(side="left", fill="x", expand=True)
        self.entry.bind("<Return>", lambda e: self.add())
        ghost_button(add, "Add", self.add, width=56).pack(side="left", padx=(6, 0))
        self.refresh()

    def add(self) -> None:
        name = self.entry.get().strip()
        if name:
            lbl = self.page.state.add_label(self.cat, name)
            self.entry.delete(0, "end")
            self.page.set_brush(self.cat, lbl["id"], toggle=False)
            self.page.state.save()

    def refresh(self) -> None:
        for w in self.chips.winfo_children():
            w.destroy()
        st, brush = self.page.state, self.page.canvas.brush
        self._chip(None, "Eraser", None, brush == (self.cat, "erase"))
        for lbl in st.cats[self.cat]["labels"]:
            self._chip(lbl, lbl["name"], lbl["color"], brush == (self.cat, lbl["id"]), st.count(self.cat, lbl["id"]))
        if brush and brush[0] == self.cat:
            name = "Eraser" if brush[1] == "erase" else st.label_name(self.cat, brush[1])
            self.sub.configure(text=f"painting · {name}", text_color=T.ACCENT)
        else:
            self.sub.configure(text="click a label to paint", text_color=T.MUTED)

    def _chip(self, lbl, text, color, selected, count=None) -> None:
        row = ctk.CTkFrame(self.chips, height=34, corner_radius=8, border_width=1,
                           fg_color=T.ACCENT_SOFT if selected else T.SURFACE,
                           border_color=T.ACCENT if selected else T.BORDER)
        row.pack(fill="x", padx=4, pady=2)
        sw = tk.Canvas(row, width=14, height=14, highlightthickness=0, bd=0,
                       bg=T.pick(T.ACCENT_SOFT if selected else T.SURFACE))
        sw.create_oval(1, 1, 13, 13, fill=color or T.pick(T.SURFACE), outline=T.pick(T.MUTED) if not color else "")
        sw.pack(side="left", padx=(10, 8), pady=9)
        name = ctk.CTkLabel(row, text=text, anchor="w", font=T.font(12), text_color=T.TEXT)
        name.pack(side="left", fill="x", expand=True)
        widgets = [row, sw, name]
        if lbl is not None:
            x = ctk.CTkLabel(row, text="✕", width=24, font=T.font(11), text_color=T.MUTED, cursor="hand2")
            x.pack(side="right", padx=(0, 6))
            x.bind("<Button-1>", lambda e, l=lbl: self.page.delete_label(self.cat, l["id"]))
            ctk.CTkLabel(row, text=str(count), font=T.font(11), text_color=T.MUTED).pack(side="right", padx=4)
            name.bind("<Double-Button-1>", lambda e, l=lbl: self.page.rename_label(self.cat, l))
        label_id = "erase" if lbl is None else lbl["id"]
        for w in widgets:
            w.bind("<Button-1>", lambda e, lid=label_id: self.page.set_brush(self.cat, lid))


class PlatePage(ctk.CTkFrame):
    title = "Plate map"
    subtitle = "Label each well with its cell type, treatment and antibody"

    def __init__(self, master, app):
        super().__init__(master, fg_color="transparent")
        self.app = app
        self.state = PlateState()

        toolbar = ctk.CTkFrame(self, fg_color="transparent")
        toolbar.pack(fill="x", pady=(0, 12))
        ctk.CTkLabel(toolbar, text="Plate", font=T.font(12), text_color=T.MUTED).pack(side="left")
        self.plate_menu = ctk.CTkSegmentedButton(
            toolbar, values=[str(n) for n in (24, 48, 96, 384)], command=self._set_plate, font=T.font(12),
            selected_color=T.ACCENT, selected_hover_color=T.ACCENT_HOVER, fg_color=T.SURFACE_ALT,
            unselected_color=T.SURFACE_ALT, unselected_hover_color=T.BORDER, text_color=T.TEXT)
        self.plate_menu.set(str(self.state.plate.n_wells))
        self.plate_menu.pack(side="left", padx=8)
        primary_button(toolbar, "Save CSV…", self.save_csv).pack(side="right")
        ghost_button(toolbar, "Copy CSV", self.copy_csv).pack(side="right", padx=6)
        ghost_button(toolbar, "Import CSV…", self.import_csv).pack(side="right")
        self.clear_btn = ghost_button(toolbar, "Clear all", self.clear_all)
        self.clear_btn.pack(side="right", padx=6)

        body = ctk.CTkFrame(self, fg_color="transparent")
        body.pack(fill="both", expand=True)
        left = Card(body)
        left.pack(side="left", fill="y", anchor="n")
        section_title(left, "Plate", "Click or drag to paint · click a row/column header to fill it").pack(
            fill="x", padx=18, pady=(14, 6))
        self.canvas = PlateCanvas(left, self.state, self._changed, self._hover)
        self.canvas.pack(padx=14, pady=4)
        self.hover = ctk.CTkLabel(left, text=" ", font=T.font(11), text_color=T.MUTED, anchor="w")
        self.hover.pack(fill="x", padx=18, pady=(2, 4))
        self.counts = ctk.CTkLabel(left, text="", font=T.font(11), text_color=T.MUTED, anchor="w")
        self.counts.pack(fill="x", padx=18, pady=(0, 14))

        right = ctk.CTkScrollableFrame(body, fg_color="transparent", width=360,
                                       scrollbar_button_color=T.BORDER, scrollbar_button_hover_color=T.MUTED)
        right.pack(side="left", fill="both", expand=True, padx=(16, 0))
        self.panels = {}
        for cat in CATEGORIES:
            self.panels[cat] = CategoryPanel(right, self, cat)
            self.panels[cat].pack(fill="x", pady=(0, 12))
        self._update_counts()

    # actions
    def set_brush(self, cat: str, label_id: str, toggle: bool = True) -> None:
        new = (cat, label_id)
        self.canvas.brush = None if (toggle and self.canvas.brush == new) else new
        for p in self.panels.values():
            p.refresh()

    def _changed(self, cat: str) -> None:
        self.panels[cat].refresh()
        self._update_counts()
        self.state.save()

    def _hover(self, text) -> None:
        self.hover.configure(text=text or " ")

    def _update_counts(self) -> None:
        n = self.state.plate.n_wells
        parts = [f"{len(self.state.cats[c]['byWell'])}/{n} {CATEGORY_TITLES[c].lower()}" for c in CATEGORIES]
        self.counts.configure(text="   ·   ".join(parts))

    def _set_plate(self, value: str) -> None:
        self.state.set_plate(PLATE_FORMATS[int(value)])
        self.state.save()
        self.canvas.redraw()
        self.refresh_all()

    def refresh_all(self) -> None:
        self.plate_menu.set(str(self.state.plate.n_wells))
        self.canvas.redraw()
        for p in self.panels.values():
            p.refresh()
        self._update_counts()

    def delete_label(self, cat: str, label_id: str) -> None:
        for w in self.state.delete_label(cat, label_id):
            self.canvas.paint(w)
        if self.canvas.brush == (cat, label_id):
            self.canvas.brush = None
        self._changed(cat)
        for p in self.panels.values():
            p.refresh()

    def rename_label(self, cat: str, lbl: dict) -> None:
        dlg = ctk.CTkInputDialog(text=f"Rename “{lbl['name']}” to:", title="Rename label")
        v = (dlg.get_input() or "").strip()
        if v:
            lbl["name"] = v
            self._changed(cat)

    def clear_all(self) -> None:
        if getattr(self, "_armed", False):
            self._armed = False
            self.clear_btn.configure(text="Clear all", text_color=T.TEXT)
            self.state.clear_wells()
            self.state.save()
            self.refresh_all()
            return
        self._armed = True
        self.clear_btn.configure(text="Click to confirm", text_color=T.DANGER)
        self.after(4000, lambda: (setattr(self, "_armed", False), self.clear_btn.configure(text="Clear all", text_color=T.TEXT)))

    def save_csv(self) -> None:
        path = filedialog.asksaveasfilename(title="Save plate map CSV", defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")], initialfile="plate_map.csv")
        if path:
            Path(path).write_text(self.state.csv_text(), encoding="utf-8")
            self.app.set_plate_map(path)
            self.app.toast(f"Saved {Path(path).name} - it will be used for the next run.")

    def write_csv_to(self, path: Path) -> Path:
        path.write_text(self.state.csv_text(), encoding="utf-8")
        return path

    def copy_csv(self) -> None:
        self.clipboard_clear()
        self.clipboard_append(self.state.csv_text())
        self.update()
        self.app.toast("Plate map CSV copied to the clipboard.")

    def import_csv(self) -> None:
        path = filedialog.askopenfilename(title="Import plate map CSV", filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if path:
            try:
                n = self.state.import_csv(Path(path))
            except Exception as exc:  # noqa: BLE001
                self.app.toast(f"Could not import: {exc}", kind="error")
                return
            self.state.save()
            self.refresh_all()
            self.app.set_plate_map(path)
            self.app.toast(f"Imported {n} labeled wells from {Path(path).name}.")

    def on_theme_change(self) -> None:
        self.canvas.redraw()
        for p in self.panels.values():
            p.refresh()
