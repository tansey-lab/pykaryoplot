"""Cytoband color schemas and small color utilities.

Schemas are verbatim copies from ``R/color.R``.
"""
from __future__ import annotations

from typing import Iterable

import numpy as np


CIRCOS = {
    "gneg": "#FFFFFF",
    "gpos25": "#C8C8C8",
    "gpos33": "#D2D2D2",
    "gpos50": "#C8C8C8",
    "gpos66": "#A0A0A0",
    "gpos75": "#828282",
    "gpos100": "#000000",
    "gpos": "#000000",
    "stalk": "#647FA4",
    "acen": "#D92F27",
    "gvar": "#DCDCDC",
    "border": "black",
}

ONLY_CENTROMERES = {
    "gneg": "#C8C8C8", "gpos25": "#C8C8C8", "gpos33": "#C8C8C8",
    "gpos50": "#C8C8C8", "gpos66": "#C8C8C8", "gpos75": "#C8C8C8",
    "gpos100": "#C8C8C8", "gpos": "#C8C8C8",
    "stalk": "#C8C8C8", "acen": "#D92F27", "gvar": "#C8C8C8",
    "border": "black",
}


def _grey(level: int) -> str:
    """R's grey<level> macro: 0 = black, 100 = white, linear interpolation."""
    v = int(round(255 * level / 100))
    v = max(0, min(255, v))
    return f"#{v:02X}{v:02X}{v:02X}"


def _build_biovizbase() -> dict[str, str]:
    d = {"gneg": _grey(100), "stalk": "#CD3333", "acen": "#8B2323",
         "gpos": _grey(0), "gvar": _grey(0), "border": "black"}
    # gpos1..gpos100 with grey(100-i) (white→black, matching R's biovizBase table)
    for i in range(1, 101):
        d[f"gpos{i}"] = _grey(100 - i)
    return d


BIOVIZBASE = _build_biovizbase()


SCHEMAS = {"circos": CIRCOS, "biovizbase": BIOVIZBASE, "only.centromeres": ONLY_CENTROMERES}


def get_cytoband_colors(color_table: dict | None = None,
                        color_schema: str = "circos") -> dict[str, str]:
    if color_table is not None:
        return dict(color_table)
    if color_schema not in SCHEMAS:
        raise ValueError(f"Unknown color_schema {color_schema!r}; "
                         f"available: {list(SCHEMAS)}")
    return dict(SCHEMAS[color_schema])


# ---------------------------------------------------------------------------
# Color manipulation helpers
# ---------------------------------------------------------------------------

_NAMED = {
    "black": (0, 0, 0), "white": (255, 255, 255),
    "red": (255, 0, 0), "green": (0, 255, 0), "blue": (0, 0, 255),
    "yellow": (255, 255, 0), "magenta": (255, 0, 255), "cyan": (0, 255, 255),
    "gray": (190, 190, 190), "grey": (190, 190, 190),
    "lightgray": (211, 211, 211), "lightgrey": (211, 211, 211),
    "darkgray": (169, 169, 169), "darkgrey": (169, 169, 169),
    "orange": (255, 165, 0), "purple": (160, 32, 240),
    "pink": (255, 192, 203), "brown": (165, 42, 42),
    "lightgreen": (144, 238, 144), "lightblue": (173, 216, 230),
    "transparent": (0, 0, 0, 0),
}


def to_rgba(color) -> tuple[float, float, float, float]:
    """Parse a color spec to (r,g,b,a) in 0..1."""
    if color is None:
        return (0, 0, 0, 0)
    if isinstance(color, tuple) and len(color) in (3, 4):
        if max(color) > 1:
            r, g, b = color[0] / 255, color[1] / 255, color[2] / 255
        else:
            r, g, b = color[:3]
        a = color[3] if len(color) == 4 else 1.0
        if a > 1:
            a = a / 255
        return (r, g, b, a)
    if isinstance(color, str):
        s = color.strip().lower()
        if s in _NAMED:
            rgb = _NAMED[s]
            if len(rgb) == 4:
                return (rgb[0] / 255, rgb[1] / 255, rgb[2] / 255, rgb[3] / 255)
            return (rgb[0] / 255, rgb[1] / 255, rgb[2] / 255, 1.0)
        if s.startswith("#"):
            hexs = s[1:]
            if len(hexs) == 3:
                hexs = "".join(ch * 2 for ch in hexs)
            if len(hexs) == 6:
                r = int(hexs[0:2], 16) / 255
                g = int(hexs[2:4], 16) / 255
                b = int(hexs[4:6], 16) / 255
                return (r, g, b, 1.0)
            if len(hexs) == 8:
                r = int(hexs[0:2], 16) / 255
                g = int(hexs[2:4], 16) / 255
                b = int(hexs[4:6], 16) / 255
                a = int(hexs[6:8], 16) / 255
                return (r, g, b, a)
    raise ValueError(f"Cannot parse color: {color!r}")


def darker(color, amount: float = 100) -> str:
    """Darken a color by `amount` (0..255 like R's darker)."""
    r, g, b, a = to_rgba(color)
    f = max(0.0, 1 - amount / 255)
    return _to_hex(r * f, g * f, b * f, a)


def lighter(color, amount: float = 100) -> str:
    r, g, b, a = to_rgba(color)
    f = amount / 255
    return _to_hex(r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f, a)


def transparent(color, amount: float = 0.5) -> str:
    """Apply an alpha multiplier (0=fully transparent, 1=opaque)."""
    r, g, b, a = to_rgba(color)
    return _to_hex(r, g, b, a * amount)


def _to_hex(r, g, b, a=1.0) -> str:
    R = int(round(r * 255))
    G = int(round(g * 255))
    B = int(round(b * 255))
    A = int(round(a * 255))
    if A == 255:
        return f"#{R:02X}{G:02X}{B:02X}"
    return f"#{R:02X}{G:02X}{B:02X}{A:02X}"


# ---------------------------------------------------------------------------
# Categorical color helpers
# ---------------------------------------------------------------------------

_DEFAULT_PALETTE = [
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
    "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
]


def col_by_chr(chrom, colors: Iterable[str] | dict | None = None) -> np.ndarray:
    chrom = np.asarray(chrom, dtype=object)
    uniq = list(dict.fromkeys(chrom.tolist()))  # stable order
    if isinstance(colors, dict):
        return np.array([colors.get(c, "#888888") for c in chrom], dtype=object)
    palette = list(colors) if colors is not None else _DEFAULT_PALETTE
    mapping = {c: palette[i % len(palette)] for i, c in enumerate(uniq)}
    return np.array([mapping[c] for c in chrom], dtype=object)


def col_by_value(value, colors=("#FFFFFF", "#000000"),
                 vmin: float | None = None, vmax: float | None = None) -> np.ndarray:
    v = np.asarray(value, dtype=float)
    lo = vmin if vmin is not None else float(np.nanmin(v))
    hi = vmax if vmax is not None else float(np.nanmax(v))
    if hi == lo:
        t = np.zeros_like(v)
    else:
        t = np.clip((v - lo) / (hi - lo), 0, 1)
    stops = [to_rgba(c) for c in colors]
    out = []
    n = len(stops) - 1
    for ti in t:
        pos = ti * n
        i = min(int(pos), n - 1)
        f = pos - i
        a = stops[i]
        b = stops[i + 1]
        out.append(_to_hex(a[0] + (b[0] - a[0]) * f,
                           a[1] + (b[1] - a[1]) * f,
                           a[2] + (b[2] - a[2]) * f,
                           a[3] + (b[3] - a[3]) * f))
    return np.array(out, dtype=object)
