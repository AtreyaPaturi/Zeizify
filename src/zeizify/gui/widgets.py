"""Small reusable widgets that carry the app's visual language."""

from __future__ import annotations

import customtkinter as ctk

from . import theme as T


class Card(ctk.CTkFrame):
    def __init__(self, master, **kw):
        kw.setdefault("fg_color", T.SURFACE)
        kw.setdefault("corner_radius", 14)
        kw.setdefault("border_width", 1)
        kw.setdefault("border_color", T.BORDER)
        super().__init__(master, **kw)


def primary_button(master, text, command, width=0, height=34, **kw) -> ctk.CTkButton:
    return ctk.CTkButton(master, text=text, command=command, width=width or 100, height=height, corner_radius=10,
                         fg_color=T.ACCENT, hover_color=T.ACCENT_HOVER, text_color=T.ON_ACCENT,
                         font=T.font(12, "bold"), **kw)


def ghost_button(master, text, command, width=0, height=34, **kw) -> ctk.CTkButton:
    return ctk.CTkButton(master, text=text, command=command, width=width or 90, height=height, corner_radius=10,
                         fg_color="transparent", hover_color=T.SURFACE_ALT, border_width=1, border_color=T.BORDER,
                         text_color=T.TEXT, font=T.font(12), **kw)


def section_title(master, title: str, subtitle: str = "") -> ctk.CTkFrame:
    f = ctk.CTkFrame(master, fg_color="transparent")
    ctk.CTkLabel(f, text=title, font=T.font(14, "bold"), text_color=T.TEXT, anchor="w").pack(fill="x")
    if subtitle:
        ctk.CTkLabel(f, text=subtitle, font=T.font(11), text_color=T.MUTED, anchor="w", justify="left").pack(fill="x")
    return f


class Toast(ctk.CTkFrame):
    """Transient message pinned to the bottom-right of the main area."""

    COLORS = {"info": T.TEXT, "success": T.SUCCESS, "error": T.DANGER, "warning": T.WARNING}

    def __init__(self, master):
        super().__init__(master, fg_color=T.SURFACE, corner_radius=12, border_width=1, border_color=T.BORDER)
        self.dot = ctk.CTkLabel(self, text="●", font=T.font(12), width=14)
        self.dot.pack(side="left", padx=(14, 6), pady=10)
        self.label = ctk.CTkLabel(self, text="", font=T.font(12), text_color=T.TEXT, wraplength=420, justify="left")
        self.label.pack(side="left", padx=(0, 16), pady=10)
        self._after = None

    def show(self, text: str, kind: str = "info", ms: int = 3500) -> None:
        self.dot.configure(text_color=self.COLORS.get(kind, T.TEXT))
        self.label.configure(text=text)
        self.place(relx=1.0, rely=1.0, x=-24, y=-24, anchor="se")
        self.lift()
        if self._after:
            self.after_cancel(self._after)
        self._after = self.after(ms, self.place_forget)
