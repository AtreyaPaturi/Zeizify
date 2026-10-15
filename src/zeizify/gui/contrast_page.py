"""
Contrast page - pick, per antibody panel, the reference field whose min/max becomes the
display range for every other field of that panel (e.g. the mock-treated well), optionally
type a custom range for any channel, and choose the figure layout: which single-channel,
nuclear-merge and co-localization panels to include and each channel's pseudocolor.
Writes intensity_report.csv and the layout PNGs.
"""

from __future__ import annotations

import threading
import time
import tkinter as tk
from pathlib import Path
from tkinter import ttk

import customtkinter as ctk

from .. import pipeline
from ..contrast import default_control_choice, fields_by_antibody
from ..layout import COLORS, candidate_panels, default_color_name, default_panels, nuclear_channel, stains_by_channel
from ..plate import parse_well
from . import theme as T
from .widgets import Card, ghost_button, primary_button, section_title


def _well_key(w: str):
    r, c = parse_well(w)
    return (r, c)


class AntibodyCard(Card):
    def __init__(self, master, antibody: str, fields: dict, channels: list[str], default: tuple[str, str], color: str,
                 stains: "dict[str, str] | None" = None):
        super().__init__(master)
        self.antibody, self.channels = antibody, channels
        head = ctk.CTkFrame(self, fg_color="transparent")
        head.pack(fill="x", padx=18, pady=(14, 8))
        sw = tk.Canvas(head, width=12, height=12, highlightthickness=0, bd=0, bg=T.pick(T.SURFACE))
        sw.create_oval(1, 1, 11, 11, fill=color, outline="")
        sw.pack(side="left", padx=(0, 8))
        ctk.CTkLabel(head, text=antibody or "(no antibody)", font=T.font(14, "bold"), text_color=T.TEXT).pack(side="left")
        stains = stains or {}
        if len(stains) > 1:
            ctk.CTkLabel(head, text="co-localization", font=T.font(10, "bold"), text_color=T.ACCENT,
                         fg_color=T.ACCENT_SOFT, corner_radius=8, height=20).pack(side="left", padx=8, ipadx=4)
        self.picked = ctk.CTkLabel(head, text="", font=T.font(11), text_color=T.MUTED)
        self.picked.pack(side="right")

        self.entries = sorted(fields.values(), key=lambda e: (_well_key(e["well"]), e["field"]))
        cols = ["field", "treatment"] + channels
        self.tree = ttk.Treeview(self, columns=cols, show="headings", style="Zz.Treeview", selectmode="browse",
                                 height=min(6, len(self.entries)))
        self.tree.heading("field", text="WELL · FIELD")
        self.tree.heading("treatment", text="TREATMENT")
        self.tree.column("field", width=120, anchor="w")
        self.tree.column("treatment", width=130, anchor="w")
        for c in channels:
            self.tree.heading(c, text=f"{c.upper()} · {stains[c]}" if stains.get(c) else c.upper())
            self.tree.column(c, width=120, anchor="center")
        self.iid_key = {}
        for i, e in enumerate(self.entries):
            iid = f"r{i}"
            self.iid_key[iid] = (e["well"], e["field"])
            vals = [f"{e['well']} · {e['field']}", e["treatment"]]
            vals += [f"{lo} – {hi}" if (lo, hi) != (None, None) else "" for lo, hi in
                     (e["channels"].get(c, (None, None)) for c in channels)]
            self.tree.insert("", "end", iid=iid, values=vals)
        self.tree.pack(fill="x", padx=18)
        default_iid = next((i for i, k in self.iid_key.items() if k == default), "r0")
        self.tree.selection_set(default_iid)
        self.tree.see(default_iid)
        self.tree.bind("<<TreeviewSelect>>", lambda e: self._update_placeholders())

        custom = ctk.CTkFrame(self, fg_color="transparent")
        custom.pack(fill="x", padx=18, pady=(10, 14))
        ctk.CTkLabel(custom, text="Custom range", font=T.font(11, "bold"), text_color=T.MUTED).pack(side="left", padx=(0, 10))
        self.custom: dict[str, tuple[ctk.CTkEntry, ctk.CTkEntry]] = {}
        for c in channels:
            box = ctk.CTkFrame(custom, fg_color=T.SURFACE_ALT, corner_radius=8)
            box.pack(side="left", padx=(0, 8))
            ctk.CTkLabel(box, text=c, font=T.font(11), text_color=T.MUTED).pack(side="left", padx=(8, 4))
            lo = ctk.CTkEntry(box, width=62, height=26, corner_radius=6, border_width=0, fg_color=T.SURFACE, font=T.font(11))
            hi = ctk.CTkEntry(box, width=62, height=26, corner_radius=6, border_width=0, fg_color=T.SURFACE, font=T.font(11))
            lo.pack(side="left", pady=4)
            ctk.CTkLabel(box, text="–", text_color=T.MUTED).pack(side="left", padx=2)
            hi.pack(side="left", padx=(0, 6), pady=4)
            self.custom[c] = (lo, hi)
        self._update_placeholders()

    def selected(self) -> tuple[str, str]:
        sel = self.tree.selection()
        return self.iid_key[sel[0]] if sel else self.iid_key["r0"]

    def _update_placeholders(self) -> None:
        well, field = self.selected()
        self.picked.configure(text=f"reference · {well} {field}")
        e = next(x for x in self.entries if (x["well"], x["field"]) == (well, field))
        for c, (lo, hi) in self.custom.items():
            a, b = e["channels"].get(c, ("", ""))
            lo.configure(placeholder_text=str(a))
            hi.configure(placeholder_text=str(b))

    def custom_ranges(self) -> tuple[dict, list[str]]:
        out, problems = {}, []
        for c, (lo_e, hi_e) in self.custom.items():
            lo_t, hi_t = lo_e.get().strip(), hi_e.get().strip()
            if not lo_t and not hi_t:
                continue
            try:
                lo, hi = int(lo_t), int(hi_t)
                if lo >= hi:
                    raise ValueError
                out[(self.antibody, c)] = (lo, hi)
            except ValueError:
                problems.append(f"{self.antibody} {c}: need integer min < max")
        return out, problems


class FigureLayoutCard(Card):
    """Choose which panels go into each field's layout PNG, and every channel's pseudocolor."""

    GROUPS = [("Single channels", lambda p, nuc: len(p) == 1),
              ("With nuclear stain", lambda p, nuc: len(p) == 2 and nuc in p),
              ("Co-localization", lambda p, nuc: len([c for c in p if c != nuc]) >= 2)]

    def __init__(self, master, channels: list[str], on_change):
        super().__init__(master)
        self.channels, self.on_change = channels, on_change
        self.nuc = nuclear_channel(channels)
        section_title(self, "Figure layout",
                      "Tick the panels to include in each field's layout PNG, left to right. Merges add the channels' "
                      "colors, so co-localized signal appears as the mixed color (e.g. green + magenta = white).").pack(
            fill="x", padx=18, pady=(14, 8))

        colors_row = ctk.CTkFrame(self, fg_color="transparent")
        colors_row.pack(fill="x", padx=18, pady=(0, 6))
        ctk.CTkLabel(colors_row, text="Colors", font=T.font(11, "bold"), text_color=T.MUTED).pack(side="left", padx=(0, 10))
        self.color_menus: dict[str, ctk.CTkOptionMenu] = {}
        for c in channels:
            box = ctk.CTkFrame(colors_row, fg_color=T.SURFACE_ALT, corner_radius=8)
            box.pack(side="left", padx=(0, 8))
            ctk.CTkLabel(box, text=c, font=T.font(11), text_color=T.TEXT).pack(side="left", padx=(10, 6))
            m = ctk.CTkOptionMenu(box, values=list(COLORS), width=96, height=26, corner_radius=6, font=T.font(11),
                                  fg_color=T.SURFACE, button_color=T.SURFACE, button_hover_color=T.BORDER,
                                  text_color=T.TEXT, dropdown_font=T.font(11), command=lambda v: self._changed())
            m.set(default_color_name(c))
            m.pack(side="left", padx=(0, 4), pady=4)
            self.color_menus[c] = m
        presets = ctk.CTkFrame(self, fg_color="transparent")
        presets.pack(fill="x", padx=18, pady=(0, 8))
        ctk.CTkLabel(presets, text="Presets", font=T.font(11, "bold"), text_color=T.MUTED).pack(side="left", padx=(0, 10))
        ghost_button(presets, "Green / Magenta (colorblind-safe)", lambda: self._preset(["Green", "Magenta", "Cyan", "Yellow"]),
                     width=230, height=28).pack(side="left", padx=(0, 6))
        ghost_button(presets, "Green / Red", lambda: self._preset(["Green", "Red", "Magenta", "Cyan"]),
                     width=110, height=28).pack(side="left", padx=(0, 6))
        ghost_button(presets, "Defaults", self._reset, width=90, height=28).pack(side="left")

        defaults = set(default_panels(channels))
        self.checks: dict[tuple, ctk.CTkCheckBox] = {}
        cands = candidate_panels(channels)
        for title, test in self.GROUPS:
            group = [p for p in cands if test(p, self.nuc) and p not in self.checks]
            if not group:
                continue
            ctk.CTkLabel(self, text=title.upper(), font=T.font(10, "bold"), text_color=T.ACCENT, anchor="w").pack(
                fill="x", padx=18, pady=(8, 2))
            row = ctk.CTkFrame(self, fg_color="transparent")
            row.pack(fill="x", padx=14)
            for i, p in enumerate(group):
                cb = ctk.CTkCheckBox(row, text=" + ".join(p), font=T.font(12), text_color=T.TEXT, checkbox_width=18,
                                     checkbox_height=18, corner_radius=5, fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER,
                                     border_color=T.BORDER, command=self._changed)
                cb.grid(row=i // 4, column=i % 4, sticky="w", padx=4, pady=3)
                if p in defaults:
                    cb.select()
                self.checks[p] = cb
        self.preview = ctk.CTkLabel(self, text="", font=T.font(11), text_color=T.MUTED, anchor="w", justify="left",
                                    wraplength=900)
        self.preview.pack(fill="x", padx=18, pady=(10, 14))
        self._update_preview()

    def _signal(self) -> list[str]:
        return [c for c in self.channels if c != self.nuc]

    def _preset(self, seq: list[str]) -> None:
        for c, color in zip(self._signal(), seq):
            self.color_menus[c].set(color)
        if self.nuc:
            self.color_menus[self.nuc].set("Blue")
        self._changed()

    def _reset(self) -> None:
        for c, m in self.color_menus.items():
            m.set(default_color_name(c))
        defaults = set(default_panels(self.channels))
        for p, cb in self.checks.items():
            cb.select() if p in defaults else cb.deselect()
        self._changed()

    def _changed(self) -> None:
        self._update_preview()
        self.on_change()

    def _update_preview(self) -> None:
        panels = self.panels()
        self.preview.configure(text=("Layout:  " + "  |  ".join("+".join(p) for p in panels)) if panels
                               else "No panels selected - layout PNGs will be skipped.",
                               text_color=T.MUTED if panels else T.WARNING)

    def panels(self) -> list[tuple[str, ...]]:
        return [p for p, cb in self.checks.items() if cb.get()]

    def colors(self) -> dict[str, str]:
        return {c: m.get() for c, m in self.color_menus.items()}


def _fmt_secs(sec: float) -> str:
    sec = int(round(sec))
    return f"{sec // 60}:{sec % 60:02d}" if sec >= 60 else f"{sec}s"


class ExportPanel(Card):
    """Live view of Apply & export: step, progress bar, current field, counts, time left, and a log."""

    STEPS = ["Intensity report", "Layout PNGs", "Run manifest"]

    def __init__(self, master, on_stop):
        super().__init__(master)
        head = ctk.CTkFrame(self, fg_color="transparent")
        head.pack(fill="x", padx=18, pady=(14, 4))
        self.title = ctk.CTkLabel(head, text="Exporting…", font=T.font(14, "bold"), text_color=T.TEXT)
        self.title.pack(side="left")
        self.step = ctk.CTkLabel(head, text="", font=T.font(10, "bold"), text_color=T.ACCENT, fg_color=T.ACCENT_SOFT,
                                 corner_radius=8, height=20)
        self.step.pack(side="left", padx=10, ipadx=4)
        self.stop_btn = ghost_button(head, "Stop", on_stop, width=70, height=28)
        self.stop_btn.pack(side="right")
        self.timing = ctk.CTkLabel(head, text="", font=T.font(11), text_color=T.MUTED)
        self.timing.pack(side="right", padx=10)

        row = ctk.CTkFrame(self, fg_color="transparent")
        row.pack(fill="x", padx=18, pady=(4, 2))
        self.bar = ctk.CTkProgressBar(row, mode="determinate", progress_color=T.ACCENT, fg_color=T.SURFACE_ALT, height=8)
        self.bar.pack(side="left", fill="x", expand=True)
        self.bar.set(0)
        self.pct = ctk.CTkLabel(row, text="0 %", width=48, font=T.font(12, "bold"), text_color=T.TEXT, anchor="e")
        self.pct.pack(side="left", padx=(10, 0))

        self.current = ctk.CTkLabel(self, text="", font=T.font(12), text_color=T.TEXT, anchor="w")
        self.current.pack(fill="x", padx=18)
        self.counts = ctk.CTkLabel(self, text="", font=T.font(11), text_color=T.MUTED, anchor="w")
        self.counts.pack(fill="x", padx=18, pady=(0, 6))
        self.log = ctk.CTkTextbox(self, height=130, font=T.mono(11), fg_color=T.SURFACE_ALT, border_width=0,
                                  corner_radius=10, text_color=T.TEXT, wrap="none")
        self.log.pack(fill="x", padx=18, pady=(0, 14))
        self.log.configure(state="disabled")
        self.n = {"written": 0, "kept": 0, "skipped": 0}
        self.t0 = self.t_layout = time.monotonic()
        self._tick()

    def set_step(self, i: int) -> None:
        self.step.configure(text=f"Step {i + 1} of {len(self.STEPS)} · {self.STEPS[i]}")
        if i == 1:
            self.t_layout = time.monotonic()
        if i != 1:
            self.bar.configure(mode="indeterminate")
            self.bar.start()
        else:
            self.bar.stop()
            self.bar.configure(mode="determinate")
            self.bar.set(0)

    def add_log(self, line: str) -> None:
        for key in self.n:
            if line.startswith(key if key != "written" else "wrote"):
                self.n[key] += 1
        self.log.configure(state="normal")
        self.log.insert("end", line + "\n")
        self.log.see("end")
        self.log.configure(state="disabled")
        if line.startswith("Writing intensity"):
            self.set_step(0)
        elif line.startswith("Writing layout"):
            self.set_step(1)

    def progress(self, done: int, total: int, name: str) -> None:
        frac = done / max(total, 1)
        self.bar.set(frac)
        self.pct.configure(text=f"{frac * 100:.0f} %")
        self.current.configure(text=f"{name}   ({done} / {total} fields)")
        self.counts.configure(text=f"{self.n['written']} written  ·  {self.n['kept']} kept (already existed)  ·  "
                                   f"{self.n['skipped']} skipped")
        el = time.monotonic() - self.t_layout
        self._eta = el / done * (total - done) if done else None

    def _tick(self) -> None:
        if not self.winfo_exists() or getattr(self, "_finished", False):
            return
        el = time.monotonic() - self.t0
        eta = getattr(self, "_eta", None)
        self.timing.configure(text=f"elapsed {_fmt_secs(el)}" + (f"  ·  about {_fmt_secs(eta)} left" if eta else ""))
        self.after(500, self._tick)

    def finish(self, title: str, kind: str) -> None:
        self._finished = True
        self.bar.stop()
        self.bar.configure(mode="determinate")
        if kind == "success":
            self.bar.set(1)
            self.pct.configure(text="100 %")
        self.title.configure(text=title, text_color={"success": T.SUCCESS, "warning": T.WARNING, "error": T.DANGER}[kind])
        self.step.pack_forget()
        self.stop_btn.pack_forget()
        self.timing.configure(text=f"took {_fmt_secs(time.monotonic() - self.t0)}")


class ContrastPage(ctk.CTkFrame):
    title = "Contrast"
    subtitle = "Pick a reference field per antibody so wells are visually comparable"

    def __init__(self, master, app):
        super().__init__(master, fg_color="transparent")
        self.app = app
        self.result = None
        self.out: "Path | None" = None
        self.body = ctk.CTkScrollableFrame(self, fg_color="transparent", scrollbar_button_color=T.BORDER,
                                           scrollbar_button_hover_color=T.MUTED)
        self.body.pack(fill="both", expand=True)
        self._empty()

    def _empty(self) -> None:
        f = ctk.CTkFrame(self.body, fg_color="transparent")
        f.pack(fill="x", pady=(140, 0))
        ctk.CTkLabel(f, text="◐", font=T.font(30), text_color=T.ACCENT).pack()
        ctk.CTkLabel(f, text="Nothing to adjust yet", font=T.font(20, "bold"), text_color=T.TEXT).pack(pady=(6, 4))
        ctk.CTkLabel(f, text="Finish an extraction on the Run page, then choose reference ranges here.",
                     font=T.font(13), text_color=T.MUTED).pack()

    def load(self, result: "pipeline.RunResult", out: Path) -> None:
        self.result, self.out = result, out
        for w in self.body.winfo_children():
            w.destroy()
        T.style_ttk(self)
        stats = result.stats
        channels = list(dict.fromkeys(s.channel_name for s in sorted(stats, key=lambda s: s.channel)))
        groups = fields_by_antibody(stats)
        defaults = default_control_choice(stats)
        antibodies = sorted(groups)
        pal = T.CATEGORY_PALETTES["antibody"] + T.CATEGORY_PALETTES["cellType"] + T.CATEGORY_PALETTES["treatment"]
        colors = {ab: pal[i % len(pal)] for i, ab in enumerate(antibodies)}

        overview = Card(self.body)
        overview.pack(fill="x", pady=(0, 12))
        section_title(overview, "Plate overview", "Wells colored by antibody panel").pack(fill="x", padx=18, pady=(14, 6))
        self._plate_overview(overview, {s.well: s.antibody for s in stats}, colors)
        legend = ctk.CTkFrame(overview, fg_color="transparent")
        legend.pack(fill="x", padx=18, pady=(4, 14))
        for ab in antibodies:
            ctk.CTkLabel(legend, text="●", text_color=colors[ab], font=T.font(14)).pack(side="left")
            ctk.CTkLabel(legend, text=ab or "(none)", font=T.font(12), text_color=T.TEXT).pack(side="left", padx=(2, 14))

        stains = stains_by_channel(stats)
        self.cards = []
        for ab in antibodies:
            card = AntibodyCard(self.body, ab, groups[ab], channels, defaults[ab], colors[ab], stains.get(ab))
            card.pack(fill="x", pady=(0, 12))
            self.cards.append(card)

        self.layout_card = FigureLayoutCard(self.body, channels, on_change=lambda: self.overwrite.select())
        self.layout_card.pack(fill="x", pady=(0, 12))

        self.export_panel: "ExportPanel | None" = None
        self.export_slot = ctk.CTkFrame(self.body, fg_color="transparent", height=1)
        self.export_slot.pack(fill="x")
        bar = Card(self.body)
        bar.pack(fill="x", pady=(0, 8))
        self.status = ctk.CTkLabel(bar, text="Empty custom boxes follow the reference field. Raw TIFFs are never modified.",
                                   font=T.font(11), text_color=T.MUTED, anchor="w", justify="left", wraplength=640)
        self.status.pack(side="left", padx=18, pady=14, fill="x", expand=True)
        self.apply_btn = primary_button(bar, "Apply & export", self.apply, width=140)
        self.apply_btn.pack(side="right", padx=(6, 14), pady=12)
        self.overwrite = ctk.CTkSwitch(bar, text="Regenerate existing PNGs", font=T.font(12), text_color=T.TEXT,
                                       progress_color=T.ACCENT)
        self.overwrite.pack(side="right", padx=6)

    def _plate_overview(self, master, well_ab: dict, colors: dict) -> None:
        plate = self.app.pages["run"].plate if hasattr(self.app.pages["run"], "plate") else self.app.pages["plate"].state.plate
        cell = max(16, min(32, 520 // plate.n_cols))
        cv = tk.Canvas(master, width=cell * (plate.n_cols + 1), height=cell * (plate.n_rows + 1), highlightthickness=0,
                       bg=T.pick(T.SURFACE))
        cv.pack(anchor="w", padx=18)
        small = (T.FAMILY, 8)
        for c in range(1, plate.n_cols + 1):
            cv.create_text(c * cell + cell / 2, cell / 2, text=str(c), fill=T.pick(T.MUTED), font=small)
        for ri, r in enumerate(plate.row_letters):
            cv.create_text(cell / 2, (ri + 1) * cell + cell / 2, text=r, fill=T.pick(T.MUTED), font=small)
            for c in range(1, plate.n_cols + 1):
                x0, y0 = c * cell + 3, (ri + 1) * cell + 3
                fill = colors.get(well_ab.get(f"{r}{c}"), T.pick(T.WELL_EMPTY))
                cv.create_oval(x0, y0, x0 + cell - 6, y0 + cell - 6, fill=fill, outline=T.pick(T.WELL_LINE))

    def apply(self) -> None:
        control, custom, problems = {}, {}, []
        for card in self.cards:
            control[card.antibody] = card.selected()
            c, p = card.custom_ranges()
            custom.update(c)
            problems += p
        if problems:
            self.status.configure(text="Fix before continuing: " + "; ".join(problems), text_color=T.DANGER)
            return
        run = self.app.pages["run"]
        overwrite = bool(self.overwrite.get())
        panels, colors = self.layout_card.panels(), self.layout_card.colors()
        cancel = threading.Event()

        if self.export_panel is not None:
            self.export_panel.destroy()
        panel = self.export_panel = ExportPanel(self.export_slot, on_stop=lambda: (
            cancel.set(), panel.stop_btn.configure(text="Stopping…", state="disabled")))
        panel.pack(fill="x", pady=(0, 12))
        panel.set_step(0)
        self.apply_btn.configure(state="disabled", text="Exporting…", fg_color=T.BORDER, text_color=T.MUTED)
        self.status.configure(text="Exporting - progress is shown above. You can keep using the rest of the app.",
                              text_color=T.MUTED)
        self.after(50, lambda: self.body._parent_canvas.yview_moveto(1.0))

        def log(line: str) -> None:
            self.app.post(panel.add_log, line)

        def work():
            res = pipeline.finalize(self.result, self.out, control, custom, overwrite_png=overwrite, log=log,
                                    panels=panels or [()], colors=colors,
                                    progress=lambda d, t, name: self.app.post(panel.progress, d, t, name),
                                    cancel=cancel)
            self.app.post(panel.set_step, 2)
            pipeline.write_manifest(self.out, Path(run.czi), Path(run.plate_map), run.plate, run.mapping,
                                    run.flip_rows, run.flip_cols, self.result, control, custom, panels, colors)
            log("Wrote run_manifest.json.")
            return res

        def done(res):
            report, n = res
            self.apply_btn.configure(state="normal", text="Apply & export", fg_color=T.ACCENT, text_color=T.ON_ACCENT)
            stopped = cancel.is_set()
            panel.finish("Export stopped" if stopped else "Export complete", "warning" if stopped else "success")
            self.status.configure(text=f"Wrote {report.name} and {n} layout PNG(s)" + (" before stopping." if stopped else "."),
                                  text_color=T.WARNING if stopped else T.SUCCESS)
            m = run.message("Export stopped" if stopped else "Report and layouts exported",
                            kind="warning" if stopped else "success")
            m.text(f"{report}\n{n} new layout PNG(s) under each panel's png/ folder.")
            from .run_page import open_path

            acts = m.actions()
            ghost_button(acts, "Open output folder", lambda: open_path(self.out), width=150).pack(side="left")
            primary_button(acts, "Analyze nuclei & colocalization  ➜",
                           lambda: (self.app.pages["analyze"].load(self.out), self.app.show("analyze")),
                           width=250).pack(side="left", padx=6)
            self.app.toast("Export stopped." if stopped else "Export complete.", kind="warning" if stopped else "success")

        def fail(exc):
            self.apply_btn.configure(state="normal", text="Apply & export", fg_color=T.ACCENT, text_color=T.ON_ACCENT)
            panel.finish("Export failed", "error")
            panel.add_log(f"ERROR: {exc}")
            self.status.configure(text=f"Export failed: {exc}", text_color=T.DANGER)

        self.app.run_bg(work, done, fail)

    def on_theme_change(self) -> None:
        T.style_ttk(self)  # keeps the user's picks; canvases refresh on the next load
