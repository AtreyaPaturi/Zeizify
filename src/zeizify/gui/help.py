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
        "intro": "Describe what is in every well: cell type, treatment, and which antibody is in each fluorophore "
                 "channel. The extraction step uses this map to name files, group wells by staining panel and "
                 "skip empty wells.",
        "sections": [
            ("Quick start", [
                "Pick the plate format (24, 48, 96 or 384 wells) at the top.",
                "Under Fluorophores, list the channels you stained (e.g. AF488, AF594). Each one gets its own "
                "Antibody card.",
                "Type a label in any card and press Enter (or Add).",
                "The new label is selected automatically - click or drag across wells to paint it.",
                "Repeat for every category. Each well shows one wedge per category, clockwise from the top "
                "(the order is written under the plate).",
            ]),
            ("Fluorophores and co-localization", [
                "Add one fluorophore per stained channel; use Quick add for common dyes or type your own.",
                "For co-localization, paint a different antibody in each channel of the same wells - e.g. CD7 in "
                "AF488 and HLA in AF594. The well's panel becomes CD7+HLA.",
                "A well can use any subset of channels; leave a channel empty where nothing was stained.",
                "Wells are grouped by panel for output folders and contrast, so only identically stained wells "
                "are compared.",
                "Fluorophore names should match the CZI channel names (AF488, AF594…) so each image is labeled "
                "with its antibody. Matching by wavelength (e.g. 'Alexa 488') also works.",
                "Double-click a fluorophore to rename it; click ✕ to remove it (twice if wells are painted).",
                "Nuclear stains such as DAPI need no fluorophore entry.",
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
                "Hover over a well to see all of its values at the bottom of the plate.",
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
                "If the file has no well IDs, they were inferred from positions. Partially imaged plates are "
                "fine: set 'First imaged well' to the imaged well closest to A1 (e.g. A2 if you started in column "
                "2 - it is pre-filled from your plate map) and press Re-map. Use Flip rows / Flip columns until "
                "the labels match the plate, then continue.",
                "If mapping fails, the message says why and shows a diagnostic plot; check the first well and the "
                "plate format (Plate map page), then retry from the same card.",
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
    "analyze": {
        "title": "Analyze",
        "intro": "Per-nucleus quantification of speckles / nuclear bodies and colocalization, run on the TIFFs "
                 "Zeizify extracted. Results go to <output folder>/analysis/.",
        "sections": [
            ("Quick start", [
                "Choose the Zeizify output folder (it is filled in automatically after a run).",
                "Tick the channels to analyse; every pair of ticked channels is tested for colocalization.",
                "Pick a field under Preview and press Preview - check that nuclei and spots are outlined correctly.",
                "Adjust settings until the preview looks right, then press Run analysis.",
            ]),
            ("Nuclei", [
                "Cellpose (deep learning, recommended) segments touching and irregular nuclei well; Classical "
                "(threshold + watershed) needs no extra install. Auto uses Cellpose when it is installed.",
                "Set the typical nucleus diameter in µm - it sets the scale for both methods and the size window "
                "(debris and clumps are rejected).",
                "Nuclei touching the image edge are excluded (incomplete).",
                "Nuclei much brighter in DAPI than the field median (condensed: mitotic / apoptotic) or dimmer "
                "(out of focus, debris), or with saturated pixels, are flagged: red in the overlay, kept in "
                "nuclei.csv with qc_pass = False, and left out of the statistics.",
                "Nuclear masks are cached: changing only spot or colocalization settings re-runs in seconds.",
            ]),
            ("Speckles / nuclear bodies", [
                "Spot diameter: the size window in µm. PML bodies are ~0.3-1 µm, nuclear speckles ~0.5-3 µm.",
                "Noise floor: a spot must rise this many robust standard deviations (MADs) above the "
                "background-subtracted nucleoplasm. Lower finds fainter spots but may split grainy nucleoplasm "
                "into false spots; check the preview.",
                "Split touching spots separates bodies that touch at an intensity dip.",
                "Threshold per channel. Per nucleus (default) adapts to each cell's staining, but it always finds "
                "the brightest bits of a nucleus - when a treatment dissolves the bodies (e.g. SUMO2/3 after "
                "ML792) it calls grainy diffuse signal 'spots'. Fixed applies one brightness floor to every well: "
                "automatic (split of all spots into dim grains vs bodies) or a value you type. Check "
                "qc/spot_brightness.png: use Fixed only if the channel shows two separated humps (separability "
                "> 0.7); keep Per nucleus for stains whose bodies persist (e.g. PML).",
                "Preview with Fixed - automatic uses the floor from the last full run.",
                "Speckles are rarely truly absent - for 'dispersed vs speckled', also use the threshold-free "
                "readouts: intensity CV in the nucleus and spot/nucleoplasm enrichment.",
            ]),
            ("Colocalization", [
                "Pixel-based, inside each nucleus: Pearson's r, Spearman's ρ, Li's ICQ, and Manders tM1/tM2 with "
                "Costes automatic thresholds (as in Fiji Coloc 2). Untresholded M1/M2 are near 1 inside a nucleus "
                "and rarely informative.",
                "Double-positive spots per nucleus: spots present in both channels (overlapping by at least one "
                "pixel) - plotted in <panel>_spots_per_nucleus.png, column n_double in colocalization.csv.",
                "Single-positive spots per nucleus: spots of one channel touching no spot of the other (n_a_only, "
                "n_b_only) - plotted with the double-positive spots in <panel>_single_double_positive.png.",
                "Object-based: fraction of spots in one channel overlapping spots in the other, overlap relative to "
                "chance (enrichment > 1 = more than random), and nearest-spot distances.",
                "Costes randomizations > 0 adds a per-nucleus significance test (slower).",
                "Wells stained with only one antibody are reported as single-stain controls - use them to check "
                "bleed-through between channels.",
            ]),
            ("Statistics", [
                "Each metric is summarised per replicate (median of its nuclei) and replicates are compared: "
                "Welch ANOVA and Games-Howell post hoc, with Benjamini-Hochberg q-values (SuperPlots, Lord 2020).",
                "A replicate is a well by default. Add a 'Replicate' column to the plate map for biological "
                "replicates; 'Field' treats fields as technical replicates and overstates certainty.",
                "With one well per condition no p-values are possible; effect sizes (Cliff's delta on nuclei) "
                "are still reported.",
                "SuperPlots show every nucleus (small dots, colored by replicate), each replicate's median (large "
                "dot, area proportional to its number of nuclei) and the mean ± SD of the replicate medians. "
                "All cell types (e.g. WT and KO) share each graph, side by side in their own color, and are also "
                "compared directly within each treatment (stats_pairwise.csv, compare = cell_type).",
                "Order: treatments run control first, then in plate-map order (e.g. MOCK, ML792, MG132, "
                "ML792+MG132); cell types wild type first. Type an order to override it.",
                "<panel>_spot_area_per_nucleus.png: average spot area per nucleus, measured at half of each spot's "
                "peak so that brighter spots do not look bigger.",
                "<panel>_spots_per_nucleus.png: spots per nucleus as a SuperPlot, plus the fraction of nuclei with "
                "0, 1-2, 3-5, ... spots (spot_count_distribution.csv).",
            ]),
        ],
    },
    "contrast": {
        "title": "Contrast",
        "intro": "Choose one reference field per antibody panel. Its brightness range is applied to every other "
                 "field of that panel, so treated wells can be compared fairly with the control. Then choose how "
                 "the figure layouts look.",
        "sections": [
            ("Quick start", [
                "Each antibody panel (e.g. CD7, or CD7+HLA for co-localization) has its own table: one row per well "
                "and field, with each channel's min – max. Column headings show which antibody is in each channel.",
                "Click the row to use as the reference (usually a control well). A control-like treatment is "
                "pre-selected for you.",
                "Press Apply & export.",
            ]),
            ("Custom ranges", [
                "Below each table you can type your own min and max for any channel.",
                "Leave both boxes empty to follow the reference field (its values are shown faintly).",
                "A custom range overrides the reference for that channel only. Min must be smaller than max.",
            ]),
            ("Figure layout and co-localization", [
                "Tick the panels each layout PNG should contain, left to right: single channels, channels merged "
                "with the nuclear stain, and co-localization merges (e.g. AF488 + AF594, with or without DAPI).",
                "Every combination of your channels is offered; the line under the checkboxes previews the order.",
                "Pick a color per channel. In merges, colors add up: overlapping green and magenta signal appears "
                "white, green and red appears yellow.",
                "Green / Magenta is the colorblind-safe choice for co-localization figures; Defaults restores "
                "the name-based colors and panels.",
                "Changing the layout turns on 'Regenerate existing PNGs' so the new figures replace the old.",
            ]),
            ("What gets written", [
                "intensity_report.csv - every field's own min/max, the antibody in each channel, the saturation "
                "flag and the range chosen here.",
                "<panel>/png/…_layout.png - the chosen panels side by side, labeled with antibody names.",
                "The raw TIFFs are never changed. Apply the ranges in Fiji or your analysis from the report.",
            ]),
            ("Tips", [
                "After Apply & export, a progress panel shows the current step (report, layout PNGs, manifest), "
                "the field being written, counts, elapsed time and an estimate of the time left, plus a live log. "
                "Stop finishes the current image and halts; images already written are kept.",
                "Existing layout PNGs are kept on re-export. Turn on 'Regenerate existing PNGs' after changing "
                "the contrast.",
                "The plate overview at the top shows which wells belong to which panel.",
                "Your choices, including panels and colors, are recorded in run_manifest.json for your methods "
                "section.",
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
