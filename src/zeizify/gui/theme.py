"""Design tokens: every color is a (light, dark) pair, as customtkinter expects."""

from __future__ import annotations

import sys
import tkinter as tk
from tkinter import ttk

import customtkinter as ctk

BG = ("#F5F6FA", "#13151B")
SIDEBAR = ("#EBEDF4", "#191C24")
SURFACE = ("#FFFFFF", "#1F232D")
SURFACE_ALT = ("#F0F2F8", "#262A36")
BORDER = ("#DCE0EB", "#2E3341")
TEXT = ("#1A1D27", "#E6E8EF")
MUTED = ("#69708A", "#9299AD")
ACCENT = ("#5650D8", "#8A85F7")
ACCENT_HOVER = ("#4741C2", "#7570E8")
ACCENT_SOFT = ("#E7E6FB", "#2B2A50")
ON_ACCENT = ("#FFFFFF", "#12121F")
SUCCESS = ("#18936A", "#41C995")
WARNING = ("#B87400", "#E5A93A")
DANGER = ("#C23B36", "#EE6A62")
WELL_EMPTY = ("#E3E6EF", "#2A2F3C")
WELL_LINE = ("#C5CBDA", "#3A4052")

FAMILY = "Segoe UI" if sys.platform == "win32" else ("SF Pro Text" if sys.platform == "darwin" else "Inter")
MONO = "Cascadia Mono" if sys.platform == "win32" else ("Menlo" if sys.platform == "darwin" else "DejaVu Sans Mono")

# Distinct, colorblind-friendlier hues per category (blues / warms / violets).
CATEGORY_PALETTES = {
    "cellType": ["#2E86AB", "#56C1D6", "#1B5E7A", "#7FB7E6", "#3D9970", "#0E9AA7", "#4E6FAE", "#8ED1C6"],
    "treatment": ["#F28E2B", "#E15759", "#F5C242", "#B5651D", "#FF9DA7", "#D4A373", "#C44536", "#EDC948"],
    "antibody": ["#8E6CEF", "#C45AB3", "#6A4C93", "#E377C2", "#B07AA1", "#5E60CE", "#D291BC", "#9D4EDD"],
}


def font(size: int = 13, weight: str = "normal") -> ctk.CTkFont:
    return ctk.CTkFont(family=FAMILY, size=size, weight=weight)


def mono(size: int = 12) -> ctk.CTkFont:
    return ctk.CTkFont(family=MONO, size=size)


def is_dark() -> bool:
    return ctk.get_appearance_mode().lower() == "dark"


def pick(pair: "tuple[str, str] | str") -> str:
    """Resolve a (light, dark) pair for raw tk widgets (Canvas, ttk) that don't do it themselves."""
    if isinstance(pair, str):
        return pair
    return pair[1] if is_dark() else pair[0]


def style_ttk(root: tk.Misc) -> None:
    """Restyle ttk.Treeview to match the current appearance mode (call again on mode change)."""
    s = ttk.Style(root)
    s.theme_use("clam")
    s.configure("Zz.Treeview", background=pick(SURFACE), fieldbackground=pick(SURFACE), foreground=pick(TEXT),
                bordercolor=pick(BORDER), lightcolor=pick(BORDER), darkcolor=pick(BORDER),
                rowheight=28, font=(FAMILY, 10), borderwidth=0)
    s.configure("Zz.Treeview.Heading", background=pick(SURFACE_ALT), foreground=pick(MUTED),
                font=(FAMILY, 9, "bold"), relief="flat", borderwidth=0, padding=(6, 6))
    s.map("Zz.Treeview", background=[("selected", pick(ACCENT_SOFT))], foreground=[("selected", pick(TEXT))])
    s.map("Zz.Treeview.Heading", background=[("active", pick(SURFACE_ALT))])
    s.layout("Zz.Treeview", [("Treeview.treearea", {"sticky": "nswe"})])
