"""
Zeizify desktop app.

Layout: a collapsible sidebar (new run, navigation, recent runs, appearance) beside a
main column with a page header and the active page. Long-running work happens on a
background thread; results are marshalled back to Tk through a queue.
"""

from __future__ import annotations

import queue
import threading
import traceback
from pathlib import Path

import customtkinter as ctk

from .. import __version__
from . import theme as T
from .contrast_page import ContrastPage
from .help import AboutDialog, HelpPanel
from .plate_page import PlatePage
from .run_page import RunPage
from .state import Settings
from .widgets import Toast

NAV = [("plate", "▦", "Plate map"), ("run", "➜", "Run"), ("contrast", "◐", "Contrast")]


class App(ctk.CTk):
    def __init__(self):
        self.settings = Settings()
        ctk.set_appearance_mode(self.settings.appearance)
        super().__init__(fg_color=T.BG)
        self.title(f"Zeizify  v{__version__}")
        self.geometry("1320x860")
        self.minsize(1080, 700)
        self._q: "queue.Queue" = queue.Queue()
        self.after(40, self._drain)
        T.style_ttk(self)

        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self._build_sidebar()

        self.main = ctk.CTkFrame(self, fg_color=T.BG, corner_radius=0)
        self.main.grid(row=0, column=1, sticky="nsew")
        header = ctk.CTkFrame(self.main, fg_color="transparent", height=64)
        header.pack(fill="x", padx=32, pady=(18, 6))
        self.expand_btn = ctk.CTkButton(header, text="☰", width=34, height=34, corner_radius=10, fg_color="transparent",
                                        hover_color=T.SURFACE_ALT, text_color=T.MUTED, font=T.font(16),
                                        command=self.toggle_sidebar)
        self.page_title = ctk.CTkLabel(header, text="", font=T.font(20, "bold"), text_color=T.TEXT, anchor="w")
        self.page_title.pack(side="left")
        self.page_sub = ctk.CTkLabel(header, text="", font=T.font(12), text_color=T.MUTED, anchor="w")
        self.page_sub.pack(side="left", padx=14, pady=(5, 0))
        self.help_btn = ctk.CTkButton(header, text="?  Help", width=84, height=32, corner_radius=16,
                                      fg_color="transparent", hover_color=T.SURFACE_ALT, border_width=1,
                                      border_color=T.BORDER, text_color=T.TEXT, font=T.font(12),
                                      command=self.toggle_help)
        self.help_btn.pack(side="right")

        content = ctk.CTkFrame(self.main, fg_color="transparent")
        content.pack(fill="both", expand=True, padx=32, pady=(6, 24))
        self.help_panel = HelpPanel(content, on_close=self.toggle_help)
        self.stack = ctk.CTkFrame(content, fg_color="transparent")
        self.stack.pack(side="left", fill="both", expand=True)
        self.current = "plate"
        self.bind("<F1>", lambda e: self.toggle_help())
        self.pages = {"plate": PlatePage(self.stack, self), "run": RunPage(self.stack, self),
                      "contrast": ContrastPage(self.stack, self)}
        self.toaster = Toast(self.main)
        self.show("run" if self.pages["plate"].state.n_labeled() else "plate")

    # ---------- sidebar ----------

    def _build_sidebar(self) -> None:
        sb = self.sidebar = ctk.CTkFrame(self, width=252, fg_color=T.SIDEBAR, corner_radius=0)
        sb.grid(row=0, column=0, sticky="nsw")
        sb.grid_propagate(False)
        sb.pack_propagate(False)

        brand = ctk.CTkFrame(sb, fg_color="transparent")
        brand.pack(fill="x", padx=16, pady=(18, 14))
        ctk.CTkLabel(brand, text="◆", font=T.font(18), text_color=T.ACCENT).pack(side="left")
        ctk.CTkLabel(brand, text="Zeizify", font=T.font(17, "bold"), text_color=T.TEXT).pack(side="left", padx=(8, 6))
        ver = ctk.CTkLabel(brand, text=f"v{__version__}", font=T.font(10, "bold"), text_color=T.ACCENT,
                           fg_color=T.ACCENT_SOFT, corner_radius=8, height=20, cursor="hand2")
        ver.pack(side="left", ipadx=4)
        ver.bind("<Button-1>", lambda e: self.show_about())
        ctk.CTkButton(brand, text="⟨", width=30, height=30, corner_radius=8, fg_color="transparent",
                      hover_color=T.BORDER, text_color=T.MUTED, font=T.font(14),
                      command=self.toggle_sidebar).pack(side="right")

        ctk.CTkButton(sb, text="＋  New run", anchor="w", height=38, corner_radius=10, fg_color=T.SURFACE,
                      hover_color=T.SURFACE_ALT, border_width=1, border_color=T.BORDER, text_color=T.TEXT,
                      font=T.font(13, "bold"), command=self.new_run).pack(fill="x", padx=12, pady=(0, 14))

        self.nav_buttons = {}
        for key, icon, label in NAV:
            b = ctk.CTkButton(sb, text=f"{icon}   {label}", anchor="w", height=36, corner_radius=10,
                              fg_color="transparent", hover_color=T.BORDER, text_color=T.TEXT, font=T.font(13),
                              command=lambda k=key: self.show(k))
            b.pack(fill="x", padx=12, pady=1)
            self.nav_buttons[key] = b

        ctk.CTkLabel(sb, text="RECENT RUNS", font=T.font(10, "bold"), text_color=T.MUTED, anchor="w").pack(
            fill="x", padx=24, pady=(22, 4))
        self.recents_box = ctk.CTkScrollableFrame(sb, fg_color="transparent", scrollbar_button_color=T.SIDEBAR,
                                                  scrollbar_button_hover_color=T.BORDER)
        self.recents_box.pack(fill="both", expand=True, padx=6)
        self.refresh_recents()

        foot = ctk.CTkFrame(sb, fg_color="transparent")
        foot.pack(fill="x", side="bottom", padx=12, pady=14)
        seg = ctk.CTkSegmentedButton(foot, values=["Light", "Dark", "System"], command=self.set_appearance,
                                     font=T.font(11), selected_color=T.ACCENT, selected_hover_color=T.ACCENT_HOVER,
                                     fg_color=T.SURFACE_ALT, unselected_color=T.SURFACE_ALT,
                                     unselected_hover_color=T.BORDER, text_color=T.TEXT)
        seg.set(self.settings.appearance.capitalize())
        seg.pack(fill="x")
        ctk.CTkButton(foot, text=f"ⓘ  About · Zeizify v{__version__}", anchor="w", height=28, corner_radius=8,
                      fg_color="transparent", hover_color=T.BORDER, text_color=T.MUTED, font=T.font(11),
                      command=self.show_about).pack(fill="x", pady=(8, 0))

    def refresh_recents(self) -> None:
        for w in self.recents_box.winfo_children():
            w.destroy()
        if not self.settings.recents:
            ctk.CTkLabel(self.recents_box, text="No runs yet", font=T.font(12), text_color=T.MUTED,
                         anchor="w").pack(fill="x", padx=12)
        for r in self.settings.recents:
            name = Path(r.get("czi", "")).name or "run"
            b = ctk.CTkButton(self.recents_box, text=name, anchor="w", height=30, corner_radius=8,
                              fg_color="transparent", hover_color=T.BORDER, text_color=T.TEXT, font=T.font(12),
                              command=lambda r=r: self.open_recent(r))
            b.pack(fill="x", padx=4, pady=1)

    def toggle_sidebar(self) -> None:
        if self.sidebar.winfo_ismapped():
            self.sidebar.grid_remove()
            self.expand_btn.pack(side="left", padx=(0, 10), before=self.page_title)
        else:
            self.sidebar.grid()
            self.expand_btn.pack_forget()

    # ---------- navigation ----------

    def show(self, key: str) -> None:
        for k, page in self.pages.items():
            page.pack_forget()
            self.nav_buttons[k].configure(fg_color=T.BORDER if k == key else "transparent",
                                          font=T.font(13, "bold" if k == key else "normal"))
        page = self.pages[key]
        page.pack(fill="both", expand=True)
        self.current = key
        if self.help_panel.winfo_ismapped():
            self.help_panel.show_page(key)
        self.page_title.configure(text=page.title)
        self.page_sub.configure(text=page.subtitle)

    def toggle_help(self) -> None:
        if self.help_panel.winfo_ismapped():
            self.help_panel.pack_forget()
            self.help_btn.configure(fg_color="transparent", text_color=T.TEXT)
        else:
            self.help_panel.show_page(self.current)
            self.help_panel.pack(side="right", fill="y", padx=(16, 0), before=self.stack)
            self.help_btn.configure(fg_color=T.ACCENT_SOFT, text_color=T.ACCENT)

    def show_about(self) -> None:
        if getattr(self, "_about", None) is not None and self._about.winfo_exists():
            self._about.focus()
            return
        self._about = AboutDialog(self)

    def new_run(self) -> None:
        self.pages["run"].reset()
        self.show("run")

    def open_recent(self, r: dict) -> None:
        self.pages["run"].load_recent(r)
        self.show("run")

    def open_contrast(self, result, out: Path) -> None:
        self.pages["contrast"].load(result, out)
        self.show("contrast")

    def set_plate_map(self, path: str) -> None:
        self.pages["run"].set_plate_map(path)

    def set_appearance(self, value: str) -> None:
        mode = value.lower()
        ctk.set_appearance_mode(mode)
        self.settings.appearance = mode
        T.style_ttk(self)
        self.after(50, lambda: [p.on_theme_change() for p in self.pages.values()])

    def toast(self, text: str, kind: str = "info") -> None:
        self.toaster.show(text, kind)

    # ---------- threading ----------

    def post(self, fn, *args) -> None:
        """Thread-safe: run fn(*args) on the Tk thread."""
        self._q.put((fn, args))

    def run_bg(self, work, on_done, on_error) -> None:
        def target():
            try:
                res = work()
            except Exception as exc:  # noqa: BLE001 - surfaced in the UI
                traceback.print_exc()
                self.post(on_error, exc)
            else:
                self.post(on_done, res)

        threading.Thread(target=target, daemon=True).start()

    def _drain(self) -> None:
        mode = ctk.get_appearance_mode()
        if mode != getattr(self, "_last_mode", mode):  # e.g. "System" and the OS theme flipped
            T.style_ttk(self)
            for p in self.pages.values():
                p.on_theme_change()
        self._last_mode = mode
        try:
            for _ in range(500):
                fn, args = self._q.get_nowait()
                fn(*args)
        except queue.Empty:
            pass
        self.after(40, self._drain)


def main() -> None:
    App().mainloop()


if __name__ == "__main__":
    main()
