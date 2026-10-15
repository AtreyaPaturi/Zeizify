"""
Run page - a conversation-style timeline of the pipeline with a composer at the bottom.

The composer collects the three inputs (CZI, plate map, output folder). Each pipeline
step posts a card to the timeline; the scene -> well QC plot arrives as a card that must
be explicitly approved (or flipped, for geometry-inferred mappings) before any well is
processed.
"""

from __future__ import annotations

import os
import subprocess
import sys
import threading
from pathlib import Path
from tkinter import filedialog

import customtkinter as ctk
from PIL import Image

from .. import pipeline
from ..mapping import MappingError
from . import theme as T
from .widgets import Card, ghost_button, primary_button


def open_path(path: Path) -> None:
    if sys.platform == "win32":
        os.startfile(path)  # noqa: S606
    else:
        subprocess.Popen(["open" if sys.platform == "darwin" else "xdg-open", str(path)])


class Message(ctk.CTkFrame):
    """Assistant-style message: small avatar on the left, content stacked to the right."""

    def __init__(self, master, title: str, kind: str = "info"):
        super().__init__(master, fg_color="transparent")
        color = {"info": T.ACCENT, "success": T.SUCCESS, "error": T.DANGER, "warning": T.WARNING}[kind]
        ctk.CTkLabel(self, text="Z", width=28, height=28, corner_radius=14, fg_color=color,
                     text_color=T.ON_ACCENT, font=T.font(12, "bold")).pack(side="left", anchor="n", padx=(0, 12), pady=2)
        self.body = ctk.CTkFrame(self, fg_color="transparent")
        self.body.pack(side="left", fill="x", expand=True)
        self.title = ctk.CTkLabel(self.body, text=title, font=T.font(13, "bold"), text_color=T.TEXT, anchor="w")
        self.title.pack(fill="x")

    def text(self, s: str, muted: bool = False) -> ctk.CTkLabel:
        lbl = ctk.CTkLabel(self.body, text=s, font=T.font(12), text_color=T.MUTED if muted else T.TEXT,
                           anchor="w", justify="left", wraplength=760)
        lbl.pack(fill="x", pady=(2, 0))
        return lbl

    def actions(self) -> ctk.CTkFrame:
        f = ctk.CTkFrame(self.body, fg_color="transparent")
        f.pack(fill="x", pady=(10, 0))
        return f


class UserBubble(ctk.CTkFrame):
    def __init__(self, master, lines: list[str]):
        super().__init__(master, fg_color="transparent")
        b = ctk.CTkFrame(self, fg_color=T.ACCENT_SOFT, corner_radius=14)
        b.pack(side="right")
        ctk.CTkLabel(b, text="\n".join(lines), font=T.font(12), text_color=T.TEXT, justify="left",
                     anchor="w").pack(padx=16, pady=10)


class PathPill(ctk.CTkButton):
    def __init__(self, master, icon: str, placeholder: str, command):
        self.icon, self.placeholder = icon, placeholder
        super().__init__(master, text=f"{icon}  {placeholder}", command=command, height=32, corner_radius=16,
                         fg_color=T.SURFACE_ALT, hover_color=T.BORDER, text_color=T.MUTED, font=T.font(12),
                         anchor="w")

    def set(self, path: "str | None") -> None:
        if path:
            self.configure(text=f"{self.icon}  {Path(path).name or path}", text_color=T.TEXT)
        else:
            self.configure(text=f"{self.icon}  {self.placeholder}", text_color=T.MUTED)


class RunPage(ctk.CTkFrame):
    title = "Run"
    subtitle = "Map scenes to wells, confirm, and extract per-well TIFFs"

    def __init__(self, master, app):
        super().__init__(master, fg_color="transparent")
        self.app = app
        self.czi: "str | None" = None
        self.plate_map: "str | None" = None
        self.out: "str | None" = None
        self.flip_rows = self.flip_cols = False
        self.first_well = "A1"
        self.mapping = None
        self.cancel = threading.Event()
        self.busy = False

        self.timeline = ctk.CTkScrollableFrame(self, fg_color="transparent", scrollbar_button_color=T.BORDER,
                                               scrollbar_button_hover_color=T.MUTED)
        self.timeline.pack(fill="both", expand=True)
        self._empty_state()

        composer = Card(self, corner_radius=18)
        composer.pack(fill="x", pady=(12, 0))
        self.hint = ctk.CTkLabel(composer, text="Choose a CZI file, a plate map and an output folder, then run.",
                                 font=T.font(13), text_color=T.MUTED, anchor="w")
        self.hint.pack(fill="x", padx=18, pady=(14, 8))
        row = ctk.CTkFrame(composer, fg_color="transparent")
        row.pack(fill="x", padx=14, pady=(0, 14))
        self.czi_pill = PathPill(row, "◧", "CZI file", self.pick_czi)
        self.czi_pill.pack(side="left")
        self.map_pill = PathPill(row, "▦", "Plate map (or use editor)", self.pick_map)
        self.map_pill.pack(side="left", padx=6)
        self.out_pill = PathPill(row, "▢", "Output folder", self.pick_out)
        self.out_pill.pack(side="left")
        self.run_btn = primary_button(row, "Run  ➜", self.run, width=96, height=36)
        self.run_btn.configure(corner_radius=18)
        self.run_btn.pack(side="right")

    # ---------- timeline helpers ----------

    def _empty_state(self) -> None:
        self.empty = ctk.CTkFrame(self.timeline, fg_color="transparent")
        self.empty.pack(fill="x", pady=(120, 0))
        ctk.CTkLabel(self.empty, text="✦", font=T.font(30), text_color=T.ACCENT).pack()
        ctk.CTkLabel(self.empty, text="What are we extracting today?", font=T.font(24, "bold"),
                     text_color=T.TEXT).pack(pady=(6, 4))
        ctk.CTkLabel(self.empty, text="Zeizify maps every scene of a plate CZI to its well, max-projects each field,\n"
                                      "and writes analysis-ready TIFFs grouped by antibody.",
                     font=T.font(13), text_color=T.MUTED, justify="center").pack()

    def _add(self, widget) -> None:
        if getattr(self, "empty", None) is not None:
            self.empty.destroy()
            self.empty = None
        widget.pack(fill="x", padx=24, pady=10)
        self.after(30, lambda: self.timeline._parent_canvas.yview_moveto(1.0))

    def message(self, title: str, kind: str = "info") -> Message:
        m = Message(self.timeline, title, kind)
        self._add(m)
        return m

    def reset(self) -> None:
        if self.busy:
            self.app.toast("A run is in progress - stop it first.", kind="warning")
            return
        for w in self.timeline.winfo_children():
            w.destroy()
        self.mapping = None
        self.flip_rows = self.flip_cols = False
        self._empty_state()

    # ---------- inputs ----------

    def pick_czi(self) -> None:
        p = filedialog.askopenfilename(title="Select the .czi file", filetypes=[("Zeiss CZI", "*.czi"), ("All files", "*.*")])
        if p:
            self.czi = p
            self.czi_pill.set(p)
            if not self.out:
                self.set_out(str(Path(p).with_name(Path(p).stem + "_zeizify")))

    def pick_map(self) -> None:
        p = filedialog.askopenfilename(title="Select the plate map CSV", filetypes=[("CSV", "*.csv"), ("All files", "*.*")])
        if p:
            self.set_plate_map(p)

    def pick_out(self) -> None:
        p = filedialog.askdirectory(title="Select an output folder")
        if p:
            self.set_out(p)

    def set_plate_map(self, p: "str | None") -> None:
        self.plate_map = p
        self.map_pill.set(p)

    def set_out(self, p: str) -> None:
        self.out = p
        self.out_pill.set(p)

    def load_recent(self, r: dict) -> None:
        self.reset()
        self.czi = r.get("czi")
        self.czi_pill.set(self.czi)
        self.set_plate_map(r.get("plate_map"))
        self.set_out(r.get("out"))

    # ---------- run flow ----------

    def _set_busy(self, busy: bool) -> None:
        self.busy = busy
        if busy:
            self.run_btn.configure(text="Stop  ■", command=self.stop, fg_color=T.DANGER, hover_color=T.DANGER)
        else:
            self.run_btn.configure(text="Run  ➜", command=self.run, fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER)

    def stop(self) -> None:
        self.cancel.set()
        self.app.toast("Stopping after the current field…", kind="warning")

    def run(self) -> None:
        if not self.czi or not Path(self.czi).exists():
            self.app.toast("Choose a CZI file first.", kind="warning")
            return
        if not self.out:
            self.app.toast("Choose an output folder first.", kind="warning")
            return
        out = Path(self.out)
        out.mkdir(parents=True, exist_ok=True)
        editor = self.app.pages["plate"]
        if not self.plate_map:
            if editor.state.n_labeled() == 0:
                self.app.toast("Label some wells in the plate map editor, or choose a plate map CSV.", kind="warning")
                return
            self.set_plate_map(str(editor.write_csv_to(out / "plate_map.csv")))
        self.plate = editor.state.plate
        self.flip_rows = self.flip_cols = False
        self.first_well = self._guess_first_well()
        self.cancel.clear()
        self.app.settings.add_recent(self.czi, self.out, self.plate_map)
        self.app.refresh_recents()

        self._add(UserBubble(self.timeline, [f"CZI  ·  {self.czi}", f"Plate map  ·  {self.plate_map}",
                                             f"Output  ·  {self.out}", f"Plate  ·  {self.plate.name}"]))
        self._map()

    def _map(self) -> None:
        m = self.message("Mapping scenes to wells…")
        bar = ctk.CTkProgressBar(m.body, mode="indeterminate", progress_color=T.ACCENT, fg_color=T.SURFACE_ALT, height=6)
        bar.pack(fill="x", pady=(8, 0))
        bar.start()
        lines: list[str] = []
        self._set_busy(True)

        def work():
            return pipeline.map_scenes(Path(self.czi), Path(self.out), self.plate, self.flip_rows, self.flip_cols,
                                       log=lines.append, first_well=self.first_well)

        def done(res):
            bar.stop()
            bar.destroy()
            self._set_busy(False)
            mapping, png = res
            self.mapping = mapping
            m.title.configure(text="Check the scene → well mapping")
            for ln in lines:
                m.text(ln, muted=True)
            self._qc_card(m, mapping, png)

        def fail(exc):
            bar.stop()
            bar.destroy()
            self._set_busy(False)
            m.title.configure(text="Mapping failed", text_color=T.DANGER)
            m.text(str(exc))
            plot = getattr(exc, "debug_plot", None) if isinstance(exc, MappingError) else None
            if plot:
                self._image(m, plot)
            if isinstance(exc, MappingError):
                m.text("If only part of the plate was imaged, enter the imaged well closest to A1 and retry. "
                       "Also check the plate format on the Plate map page.", muted=True)
                self._first_well_row(m)

        self.app.run_bg(work, done, fail)

    def _image(self, m: Message, path: Path) -> None:
        img = Image.open(path)
        w = 760
        h = int(img.height * w / img.width)
        frame = ctk.CTkFrame(m.body, fg_color=T.SURFACE, corner_radius=12, border_width=1, border_color=T.BORDER)
        frame.pack(fill="x", pady=(10, 0))
        ctk.CTkLabel(frame, text="", image=ctk.CTkImage(img, img, size=(w, h))).pack(padx=8, pady=8)
        frame.bind("<Double-Button-1>", lambda e: open_path(path))

    def _guess_first_well(self) -> str:
        """The labeled well closest to A1 in the plate map - a good default for partial plates."""
        try:
            from ..plate import load_plate_map, parse_well

            wells = [parse_well(w) for w in load_plate_map(Path(self.plate_map))]
            if wells:
                letters = self.plate.row_letters
                r = min(letters.index(w[0]) for w in wells if w[0] in letters)
                c = min(w[1] for w in wells)
                return f"{letters[r]}{c}"
        except Exception:  # noqa: BLE001 - a guess only
            pass
        return "A1"

    def _first_well_row(self, m: Message) -> None:
        row = ctk.CTkFrame(m.body, fg_color="transparent")
        row.pack(fill="x", pady=(10, 0))
        ctk.CTkLabel(row, text="First imaged well", font=T.font(12), text_color=T.TEXT).pack(side="left", padx=(0, 8))
        entry = ctk.CTkEntry(row, width=70, height=30, corner_radius=8, border_color=T.BORDER, fg_color=T.SURFACE_ALT,
                             font=T.font(12))
        entry.insert(0, self.first_well)
        entry.pack(side="left")
        ctk.CTkLabel(row, text="the imaged well closest to A1, e.g. A2 if you started in column 2", font=T.font(11),
                     text_color=T.MUTED).pack(side="left", padx=10)

        def apply():
            from ..plate import normalize_well

            try:
                self.first_well = normalize_well(entry.get())
            except ValueError:
                self.app.toast("Enter a well ID such as A2.", kind="warning")
                return
            btn.configure(state="disabled")
            self._map()

        btn = ghost_button(row, "Re-map", apply, width=80, height=30)
        btn.pack(side="left")
        entry.bind("<Return>", lambda e: apply())

    def _qc_card(self, m: Message, mapping, png: Path) -> None:
        if mapping.source == "metadata":
            m.text("Well IDs were read directly from the well names ZEN stored in the CZI, so there is no orientation "
                   "to flip. Skim the plot for anything obviously wrong, then continue.")
        else:
            row = "bottom" if self.flip_rows else "top"
            col = "right" if self.flip_cols else "left"
            m.text("This CZI has no well-ID metadata, so wells were inferred from stage positions. Confirm that the "
                   f"layout matches how the plate was scanned. Current assumption: first imaged well "
                   f"{mapping.first_well}, A1 at {row}-{col}.")
            self._first_well_row(m)
        self._image(m, png)
        a = m.actions()
        btns = []

        def choose(choice):
            for b in btns:
                b.configure(state="disabled")
            if choice == "ok":
                self._extract()
                return
            if choice == "cancel":
                self.message("Stopped before processing any wells.", kind="warning")
                return
            if choice in ("rows", "both"):
                self.flip_rows = not self.flip_rows
            if choice in ("cols", "both"):
                self.flip_cols = not self.flip_cols
            self._map()

        btns.append(primary_button(a, "Looks correct — continue", lambda: choose("ok"), width=200))
        if mapping.source != "metadata":
            for text, key in (("Flip rows", "rows"), ("Flip columns", "cols"), ("Flip both", "both")):
                btns.append(ghost_button(a, text, lambda k=key: choose(k)))
        btns.append(ghost_button(a, "Cancel", lambda: choose("cancel")))
        btns.append(ghost_button(a, "Open plot", lambda: open_path(png)))
        for b in btns:
            b.pack(side="left", padx=(0, 6))

    def _extract(self) -> None:
        m = self.message("Extracting wells…")
        status = m.text("Starting…", muted=True)
        bar = ctk.CTkProgressBar(m.body, mode="determinate", progress_color=T.ACCENT, fg_color=T.SURFACE_ALT, height=6)
        bar.set(0)
        bar.pack(fill="x", pady=(8, 6))
        log = ctk.CTkTextbox(m.body, height=180, font=T.mono(11), fg_color=T.SURFACE, border_width=1,
                             border_color=T.BORDER, corner_radius=10, text_color=T.TEXT, wrap="none")
        log.pack(fill="x")
        log.configure(state="disabled")
        self._set_busy(True)

        def append(s):
            log.configure(state="normal")
            log.insert("end", s + "\n")
            log.see("end")
            log.configure(state="disabled")

        def prog(done, total):
            bar.set(done / max(total, 1))
            status.configure(text=f"{done} / {total} fields")

        def work():
            return pipeline.extract(Path(self.czi), Path(self.plate_map), Path(self.out), self.mapping,
                                    log=lambda s: self.app.post(append, s),
                                    progress=lambda d, t: self.app.post(prog, d, t), cancel=self.cancel)

        def done(result):
            self._set_busy(False)
            self.result = result
            sat = sum(s.saturated for s in result.stats)
            ok = not result.skipped
            m.title.configure(text="Extraction finished" if not self.cancel.is_set() else "Extraction stopped")
            pipeline.write_manifest(Path(self.out), Path(self.czi), Path(self.plate_map), self.plate, self.mapping,
                                    self.flip_rows, self.flip_cols, result)
            s = self.message(f"{len(result.stats)} channel TIFFs ready", kind="success" if ok else "warning")
            s.text(f"{len({(x.well, x.field) for x in result.stats})} fields across "
                   f"{len({x.well for x in result.stats})} wells and {len({x.antibody for x in result.stats})} antibodies."
                   + (f"  {sat} channel image(s) are saturated at the sensor ceiling." if sat else ""))
            if result.skipped:
                s.text(f"{len(result.skipped)} item(s) skipped: " + "; ".join(f"{n} ({r})" for n, r in result.skipped[:6])
                       + (" …" if len(result.skipped) > 6 else ""), muted=True)
            a = s.actions()
            if result.stats:
                primary_button(a, "Choose contrast  ➜", lambda: self.app.open_contrast(result, Path(self.out)),
                               width=170).pack(side="left", padx=(0, 6))
            ghost_button(a, "Open output folder", lambda: open_path(Path(self.out)), width=150).pack(side="left")

        def fail(exc):
            self._set_busy(False)
            m.title.configure(text="Extraction failed", text_color=T.DANGER)
            m.text(str(exc))

        self.app.run_bg(work, done, fail)

    def on_theme_change(self) -> None:
        pass
