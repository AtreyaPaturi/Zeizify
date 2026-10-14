"""
Read the parts of a CZI's XML metadata that Zeizify needs - without touching pixel data.

All functions take the metadata XML *string* (``CziFile(...).metadata()``), so they are
testable without a CZI on disk. Element paths follow the ZEN ImageDocument schema and
match those used by ZEISS' own czitools (sebi06/czitools):

    Scenes   Metadata/Information/Image/Dimensions/S/Scenes/Scene   (@Index, @Name, ArrayName, Shape/RowIndex|ColumnIndex)
    Channels Metadata/Information/Image/Dimensions/Channels/Channel (@Id, @Name), in C-index order
    Scaling  Metadata/Scaling/Items/Distance[@Id='X']/Value         (meters per pixel)
    Bits     Metadata/Information/Image/ComponentBitCount           (falls back to per-channel values)
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass

_IMAGE = "Information/Image"


@dataclass(frozen=True)
class SceneInfo:
    index: int
    name: str  # field/position name assigned by ZEN, e.g. "P1"
    array_name: str  # well ID assigned by ZEN, e.g. "A1" ("" if absent)
    row_index: "int | None"  # Shape/RowIndex, if present
    col_index: "int | None"  # Shape/ColumnIndex, if present


def _metadata_root(xml: str) -> ET.Element:
    root = ET.fromstring(xml)
    if root.tag == "Metadata":
        return root
    md = root.find("Metadata")
    if md is None:
        raise ValueError("CZI XML has no <Metadata> element")
    return md


def _int_or_none(text: "str | None") -> "int | None":
    try:
        return int(text.strip()) if text is not None else None
    except ValueError:
        return None


def read_scenes(xml: str) -> list[SceneInfo]:
    """Every <Scene> under Information/Image/Dimensions/S/Scenes, keyed by its own Index
    attribute (not its position in the list - ZEN writes them in order, but the Index is
    what subblocks reference, so that is what we trust)."""
    md = _metadata_root(xml)
    out = []
    for i, sc in enumerate(md.findall(f"{_IMAGE}/Dimensions/S/Scenes/Scene")):
        idx = _int_or_none(sc.get("Index"))
        shape = sc.find("Shape")
        out.append(
            SceneInfo(
                index=i if idx is None else idx,
                name=sc.get("Name", "") or "",
                array_name=(sc.findtext("ArrayName") or "").strip(),
                row_index=_int_or_none(shape.findtext("RowIndex")) if shape is not None else None,
                col_index=_int_or_none(shape.findtext("ColumnIndex")) if shape is not None else None,
            )
        )
    return out


def read_channel_names(xml: str) -> list[str]:
    """Channel names in C-axis order. Falls back to DisplaySetting/Channels, then to the
    first-seen <Channel Name=...> anywhere, only if the Information block has none."""
    md = _metadata_root(xml)
    for path in (f"{_IMAGE}/Dimensions/Channels/Channel", "DisplaySetting/Channels/Channel"):
        chans = md.findall(path)
        if chans:
            return [c.get("Name") or c.findtext("ShortName") or f"C{i}" for i, c in enumerate(chans)]
    seen: list[str] = []
    for name in re.findall(r'<Channel\b[^>]*\bName="([^"]*)"', xml):
        if name not in seen:
            seen.append(name)
    return seen


def read_field_names(xml: str) -> dict[int, str]:
    """Scene index -> ZEN field/position name (e.g. 'P1')."""
    return {s.index: s.name for s in read_scenes(xml)}


def read_bit_ceiling(xml: str) -> "int | None":
    """Sensor ceiling, 2**bits - 1 (e.g. 16383 for 14-bit). None if absent or ambiguous."""
    md = _metadata_root(xml)
    bits = {_int_or_none(e.text) for e in md.findall(f"{_IMAGE}/ComponentBitCount")}
    if not bits - {None}:
        bits = {_int_or_none(e.text) for e in md.findall(f"{_IMAGE}/Dimensions/Channels/Channel/ComponentBitCount")}
    bits.discard(None)
    if len(bits) == 1:
        return 2 ** bits.pop() - 1
    return None


def read_pixel_size_um(xml: str) -> "float | None":
    """Physical X pixel size in micrometers. CZI stores meters per pixel."""
    md = _metadata_root(xml)
    for d in md.findall("Scaling/Items/Distance"):
        if d.get("Id") == "X":
            try:
                v = float((d.findtext("Value") or "").strip())
            except ValueError:
                return None
            return v * 1e6 if v > 0 else None
    return None
