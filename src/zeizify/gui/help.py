"""
In-app help: one guide per page, shown in a side panel (toggle with the Help button or F1),
plus the About dialog that reports the exact version a user is running.

Edit HELP below to update the guides - each page key maps to an intro and a list of
(heading, [bullet, ...]) sections.
"""

from __future__ import annotations

import platform
import webbrowser
from importlib import metadata as importlib_metadata

import customtkinter as ctk

from .. import __version__
from . import theme as T

REPO_URL = "https://github.com/AtreyaPaturi/Zeizify"

HELP: dict[str, dict] = {
    "plate": {
        "title": "Plate map",
        "intro": "Describe what is in every well. The extraction step uses this map to name files, group them "
                 "by antibody and skip empty wells.",
        "sections": [
            ("Quick start", [
                "Pick the plate format (24, 48, 96 or 384 wells) at the top.",
                "Type a label under Cell type, Treatment or Antibody and press Enter (or Add).",
                "The new label is selected automatically - click or drag across wells to paint it.",
                "Repeat for the other two categories. Each well shows all three as wedges.",
            ]),
            ("Painting", [
                "Click a well to paint it; click it again with the same label to clear it.",
                "Drag across wells to paint (or clear) a run of them in one stroke.",
                "Click a row letter, column number or 'all' to fill that whole line. Clicking again when the "
                "line is already filled clears it.",
                "Select Eraser to clear wells for that category only.",
                "Click the selected label again to stop painting.",
            ]),
            ("Managing labels", [
                "Double-click a label's name to rename it - every painted well updates.",
                "Click ✕ to delete a label; its wells become empty for that category.",
                "The number next to a label is how many wells carry it.",
                "Hover over a well to see all three values at the bottom of the plate.",
            ]),
            ("Saving and sharing", [
                "Your map autosaves - closing the app never loses work.",
                "Save CSV… writes the plate map file and selects it for the next run.",
                "Import CSV… loads an existing plate map (also hand-made ones in Excel).",
                "Copy CSV puts the map on the clipboard, e.g. to paste into a spreadsheet.",
                "Clear all empties every well (click twice to confirm); labels are kept.",
            ]),
            ("Tips", [
                "Wells with no label at all are skipped during extraction.",
                "Include the word Mock, Control, DMSO, Vehicle or Untreated in control treatments - the "
                "Contrast page then pre-selects them as the reference.",
                "Labels become part of file names; special characters are replaced by '_'.",
            ]),
        ],
    },
    "run": {
        "title": "Run",
        "intro": "Map every scene in the CZI to its well, confirm the mapping, then extract one TIFF per well, "
                 "field and channel.",
        "sections": [
            ("Quick start", [
                "In the bar at the bottom, choose the CZI file, a plate map and an output folder.",
                "Leave the plate map empty to use the map currently in the Plate map editor.",
                "Press Run.",
            ]),
            ("Checking the mapping", [
                "A plot of every field's stage position, labeled by well, appears in the timeline.",
                "If the well IDs came from the CZI metadata (most ZEN plate acquisitions), just check that the "
                "plot looks like your plate and press 'Looks correct — continue'.",
                "If the file has no well IDs, they were inferred from positions. Use Flip rows / Flip columns "
                "until A1 is where it physically was on the plate, then continue.",
                "Nothing is extracted until you confirm. Cancel stops safely.",
                "Open plot shows the full-resolution image.",
            ]),
            ("During extraction", [
                "The progress bar counts fields; the log shows each channel's min/max intensity.",
                "[SATURATED] marks images that hit the camera's maximum value - consider re-imaging or "
                "excluding them from quantification.",
                "Stop finishes the current field and halts. Run again with the same output folder to resume - "
                "finished fields are not redone.",
            ]),
            ("When it finishes", [
                "Choose contrast opens the Contrast page for this run.",
                "Open output folder shows the results: one folder per antibody containing the TIFFs.",
                "Everything about the run (inputs, settings, software versions) is saved in run_manifest.json.",
            ]),
            ("Tips", [
                "Recent runs in the sidebar refill all three inputs with one click.",
                "New run clears the timeline (inputs you picked are kept).",
                "TIFFs carry the physical pixel size, so Fiji/ImageJ shows correct µm scale bars.",
                "Very large files are fine: only one field is read into memory at a time.",
            ]),
        ],
    },
    "contrast": {
        "title": "Contrast",
        "intro": "Choose one reference field per antibody. Its brightness range is applied to every other field "
                 "of that antibody, so treated wells can be compared fairly with the control.",
        "sections": [
            ("Quick start", [
                "Each antibody has its own table: one row per well and field, with each channel's min – max.",
                "Click the row to use as the reference (usually a control well). A control-like treatment is "
                "pre-selected for you.",
                "Press Apply & export.",
            ]),
            ("Custom ranges", [
                "Below each table you can type your own min and max for any channel.",
                "Leave both boxes empty to follow the reference field (its values are shown faintly).",
                "A custom range overrides the reference for that channel only. Min must be smaller than max.",
            ]),
            ("What gets written", [
                "intensity_report.csv - every field's own min/max, saturation flag and the range chosen here.",
                "<antibody>/png/…_layout.png - each channel and its DAPI merge side by side, using these ranges.",
                "The raw TIFFs are never changed. Apply the ranges in Fiji or your analysis from the report.",
            ]),
            ("Tips", [
                "Existing layout PNGs are kept on re-export. Turn on 'Regenerate existing PNGs' after changing "
                "the contrast.",
                "The plate overview at the top shows which wells belong to which antibody.",
                "Your choices are recorded in run_manifest.json for your methods section.",
            ]),
        ],
    },
}


class HelpPanel(ctk.CTkFrame):
    """Right-hand side panel showing the guide for the current page."""

    def __init__(self, master, on_close):
        super().__init__(master, width=340, fg_color=T.SURFACE, corner_radius=14, border_width=1, border_color=T.BORDER)
        self.pack_propagate(False)
        head = ctk.CTkFrame(self, fg_color="transparent")
        head.pack(fill="x", padx=16, pady=(14, 4))
        ctk.CTkLabel(head, text="?", width=24, height=24, corner_radius=12, fg_color=T.ACCENT_SOFT,
                     text_color=T.ACCENT, font=T.font(12, "bold")).pack(side="left")
        self.title = ctk.CTkLabel(head, text="", font=T.font(14, "bold"), text_color=T.TEXT)
        self.title.pack(side="left", padx=8)
        ctk.CTkButton(head, text="✕", width=28, height=28, corner_radius=8, fg_color="transparent",
                      hover_color=T.SURFACE_ALT, text_color=T.MUTED, font=T.font(12), command=on_close).pack(side="right")
        self.body = ctk.CTkScrollableFrame(self, fg_color="transparent", scrollbar_button_color=T.BORDER,
                                           scrollbar_button_hover_color=T.MUTED)
        self.body.pack(fill="both", expand=True, padx=4, pady=(0, 8))
        foot = ctk.CTkFrame(self, fg_color="transparent")
        foot.pack(fill="x", padx=16, pady=(0, 12))
        ctk.CTkLabel(foot, text=f"Zeizify v{__version__}  ·  F1 toggles help", font=T.font(10),
                     text_color=T.MUTED).pack(side="left")
        link = ctk.CTkLabel(foot, text="Documentation ↗", font=T.font(10, "bold"), text_color=T.ACCENT, cursor="hand2")
        link.pack(side="right")
        link.bind("<Button-1>", lambda e: webbrowser.open(REPO_URL + "#readme"))

    def show_page(self, key: str) -> None:
        g = HELP[key]
        self.title.configure(text=f"Help · {g['title']}")
        for w in self.body.winfo_children():
            w.destroy()
        ctk.CTkLabel(self.body, text=g["intro"], font=T.font(12), text_color=T.TEXT, wraplength=280,
                     justify="left", anchor="w").pack(fill="x", padx=10, pady=(6, 4))
        for heading, bullets in g["sections"]:
            ctk.CTkLabel(self.body, text=heading.upper(), font=T.font(10, "bold"), text_color=T.ACCENT,
                         anchor="w").pack(fill="x", padx=10, pady=(14, 4))
            for b in bullets:
                row = ctk.CTkFrame(self.body, fg_color="transparent")
                row.pack(fill="x", padx=10, pady=2)
                ctk.CTkLabel(row, text="•", font=T.font(12), text_color=T.MUTED, width=12).pack(side="left", anchor="n")
                ctk.CTkLabel(row, text=b, font=T.font(12), text_color=T.TEXT, wraplength=262, justify="left",
                             anchor="w").pack(side="left", fill="x", expand=True)
        self.body._parent_canvas.yview_moveto(0)


def version_info() -> dict[str, str]:
    """Everything worth pasting into a bug report or a methods section."""
    info = {"Zeizify": __version__, "Python": platform.python_version(), "OS": platform.platform()}
    for pkg in ("czifile", "tifffile", "imagecodecs", "numpy", "scipy", "customtkinter"):
        try:
            info[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            pass
    return info


class AboutDialog(ctk.CTkToplevel):
    def __init__(self, master):
        super().__init__(master, fg_color=T.BG)
        self.title("About Zeizify")
        self.resizable(False, False)
        self.transient(master)
        self.after(50, self.grab_set)

        ctk.CTkLabel(self, text="◆", font=T.font(30), text_color=T.ACCENT).pack(pady=(22, 0))
        ctk.CTkLabel(self, text="Zeizify", font=T.font(22, "bold"), text_color=T.TEXT).pack()
        ctk.CTkLabel(self, text=f"Version {__version__}", font=T.font(13, "bold"), text_color=T.ACCENT,
                     fg_color=T.ACCENT_SOFT, corner_radius=10, height=26).pack(pady=(6, 4), ipadx=10)
        ctk.CTkLabel(self, text="Per-well extraction of Zeiss CZI multi-well plates", font=T.font(12),
                     text_color=T.MUTED).pack(pady=(0, 12))

        card = ctk.CTkFrame(self, fg_color=T.SURFACE, corner_radius=12, border_width=1, border_color=T.BORDER)
        card.pack(fill="x", padx=24)
        self.info = version_info()
        for k, v in self.info.items():
            if k == "Zeizify":
                continue
            row = ctk.CTkFrame(card, fg_color="transparent")
            row.pack(fill="x", padx=14, pady=1)
            ctk.CTkLabel(row, text=k, font=T.font(11), text_color=T.MUTED, anchor="w").pack(side="left")
            ctk.CTkLabel(row, text=v if len(v) < 34 else v[:31] + "…", font=T.mono(11), text_color=T.TEXT).pack(side="right")

        btns = ctk.CTkFrame(self, fg_color="transparent")
        btns.pack(padx=24, pady=(16, 22))
        self.copy_btn = ctk.CTkButton(btns, text="Copy version info", command=self.copy, height=32, corner_radius=10,
                                      fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER, text_color=T.ON_ACCENT,
                                      font=T.font(12, "bold"))
        self.copy_btn.pack(side="left", padx=4)
        ctk.CTkButton(btns, text="GitHub ↗", command=lambda: webbrowser.open(REPO_URL), height=32, width=90,
                      corner_radius=10, fg_color="transparent", hover_color=T.SURFACE_ALT, border_width=1,
                      border_color=T.BORDER, text_color=T.TEXT, font=T.font(12)).pack(side="left", padx=4)

    def copy(self) -> None:
        self.clipboard_clear()
        self.clipboard_append("\n".join(f"{k}: {v}" for k, v in self.info.items()))
        self.copy_btn.configure(text="Copied ✓")
