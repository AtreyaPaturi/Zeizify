"""Keep every place that states the version in sync with zeizify.__version__."""

import re
from pathlib import Path

import zeizify

ROOT = Path(__file__).resolve().parents[1]


def test_citation_version_matches():
    cff = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    assert re.search(r"^version:\s*(\S+)", cff, re.M).group(1) == zeizify.__version__


def test_changelog_has_entry_for_version():
    assert f"## [{zeizify.__version__}]" in (ROOT / "CHANGELOG.md").read_text(encoding="utf-8")


def test_every_page_has_help():
    from zeizify.gui.help import HELP

    assert set(HELP) == {"plate", "run", "contrast"}
    assert all(g["intro"] and g["sections"] for g in HELP.values())
