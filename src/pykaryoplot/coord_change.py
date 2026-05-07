"""Genomic-to-plot coordinate transforms.

Direct port of ``R/getCoordChangeFunctions.R``. Two layouts:

* **stacked** (plot.types 1, 2, 6): chromosomes drawn one above the other
  with a fixed per-chromosome height.
* **single-line** (plot.types 3, 4, 5, 7): all chromosomes in a single row,
  separated by ``ideogramlateralmargin``.

The returned ``coord_change`` callable takes per-element ``chrom`` plus
optional ``x`` (in bp) and ``y`` (in data-panel units) arrays and returns
plot coordinates as ``(x_plot, y_plot)`` numpy arrays of the same length.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

import numpy as np

from .plot_params import SINGLE_LINE_TYPES


def _chr_height(pp: dict) -> float:
    return (pp["data2outmargin"] + pp["data2height"] + pp["data2inmargin"]
            + pp["ideogramheight"]
            + pp["data1inmargin"] + pp["data1height"] + pp["data1outmargin"])


def _ideogram_mid_stacked(chr_names: list[str], pp: dict) -> dict[str, float]:
    """For stacked layouts: y of each chromosome's ideogram midline."""
    chr_height = _chr_height(pp)
    n = len(chr_names)
    # In R: `chrs <- c(length(chr.names):1)` — chromosome at top is index 1
    # of the reversed numbering. So chr 1 (first listed) is at the top.
    out = {}
    for i, name in enumerate(chr_names):
        chr_num = n - i  # top chromosome -> n, bottom -> 1
        out[name] = (pp["bottommargin"] + (chr_num - 1) * chr_height
                     + pp["data2outmargin"] + pp["data2height"] + pp["data2inmargin"]
                     + pp["ideogramheight"] / 2.0)
    return out


def _ideogram_mid_single_line(pp: dict) -> float:
    return (pp["bottommargin"] + pp["data2outmargin"] + pp["data2height"]
            + pp["data2inmargin"] + pp["ideogramheight"] / 2.0)


def _scale_y(y: np.ndarray, panel, pp: dict, chr_mid: np.ndarray) -> np.ndarray:
    """Scale data-panel y-values to plot y."""
    if panel == "ideogram":
        rng = pp["dataideogrammax"] - pp["dataideogrammin"]
        ys = (y - pp["dataideogrammin"]) / rng * pp["ideogramheight"]
        return chr_mid - pp["ideogramheight"] / 2 + ys
    if panel == "all":
        all_h = _chr_height(pp) - pp["data1outmargin"] - pp["data2outmargin"]
        rng = pp["dataallmax"] - pp["dataallmin"]
        ys = (y - pp["dataallmin"]) / rng * all_h
        return chr_mid - pp["ideogramheight"] / 2 - pp["data2inmargin"] - pp["data2height"] + ys
    if panel == 1:
        rng = pp["data1max"] - pp["data1min"]
        ys = (y - pp["data1min"]) / rng * pp["data1height"]
        return chr_mid + pp["ideogramheight"] / 2 + pp["data1inmargin"] + ys
    if panel == 2:
        rng = pp["data2max"] - pp["data2min"]
        ys = (y - pp["data2min"]) / rng * pp["data2height"]
        if pp.get("data2_invert", True):
            return chr_mid - pp["ideogramheight"] / 2 - pp["data2inmargin"] - ys
        # Non-inverted: panel 2 grows upward like panel 1, sitting just below
        # the ideogram. Data2max is the value closest to the ideogram.
        return (chr_mid - pp["ideogramheight"] / 2
                - pp["data2inmargin"] - pp["data2height"] + ys)
    raise ValueError(f"Invalid data.panel: {panel}")


@dataclass
class CoordChange:
    """Bundle of coordinate-mapping helpers attached to a KaryoPlot."""
    coord_change: Callable
    ideogram_mid: Callable[[str | np.ndarray], np.ndarray]
    chr_height: float
    plot_xlim: tuple[float, float]
    plot_ylim: tuple[float, float]


def build_coord_change(plot_type: int, chr_names: list[str],
                       chr_starts: dict[str, int], chr_ends: dict[str, int],
                       pp: dict) -> CoordChange:
    """Construct the coordinate-change callable for a karyoplot."""
    chr_height = _chr_height(pp)
    n_chr = len(chr_names)

    if plot_type in SINGLE_LINE_TYPES:
        return _build_single_line(plot_type, chr_names, chr_starts, chr_ends, pp, chr_height)
    return _build_stacked(plot_type, chr_names, chr_starts, chr_ends, pp, chr_height, n_chr)


def _build_stacked(plot_type, chr_names, chr_starts, chr_ends, pp, chr_height, n_chr):
    mid_map = _ideogram_mid_stacked(chr_names, pp)
    chr_lens = {c: chr_ends[c] - chr_starts[c] for c in chr_names}
    max_chr_len = max(chr_lens.values()) if chr_lens else 1
    genome_width = 1 - pp["leftmargin"] - pp["rightmargin"]

    def ideogram_mid(chrom):
        chrom = np.atleast_1d(np.asarray(chrom, dtype=object))
        return np.array([mid_map[c] for c in chrom], dtype=float)

    def coord_change(chrom=None, x=None, y=None, data_panel=1):
        # plot_type-specific data_panel coercion
        if plot_type == 1 and (data_panel is None or data_panel == 2):
            data_panel = 1
        if plot_type == 6 and data_panel in (None, 1, 2):
            data_panel = "ideogram"

        x_plot = None
        y_plot = None

        if x is not None and chrom is not None:
            x = np.asarray(x, dtype=float)
            chrom_arr = np.asarray(chrom, dtype=object)
            starts = np.array([chr_starts[c] for c in chrom_arr], dtype=float)
            x_plot = pp["leftmargin"] + ((x - starts) / max_chr_len) * genome_width

        if chrom is not None:
            chrom_arr = np.asarray(chrom, dtype=object)
            mid = np.array([mid_map[c] for c in chrom_arr], dtype=float)
            if y is None:
                y_plot = mid
            else:
                y_plot = _scale_y(np.asarray(y, dtype=float), data_panel, pp, mid)

        return x_plot, y_plot

    # plot region bounds — matches plotKaryotype.R lines 340-348
    if plot_type in (1, 2, 6):
        ylim = (0, pp["bottommargin"] + n_chr * chr_height + pp["topmargin"])
    else:
        ylim = (0, pp["bottommargin"] + chr_height + pp["topmargin"])
    return CoordChange(coord_change, ideogram_mid, chr_height, (0.0, 1.0), ylim)


def _build_single_line(plot_type, chr_names, chr_starts, chr_ends, pp, chr_height):
    n = len(chr_names)
    chr_idx = {c: i for i, c in enumerate(chr_names)}
    chr_lens = {c: chr_ends[c] - chr_starts[c] for c in chr_names}
    all_chr_len = sum(chr_lens.values()) if chr_lens else 1
    # cumulative sum of preceding chromosome lengths
    prev = {}
    acc = 0
    for c in chr_names:
        prev[c] = acc
        acc += chr_lens[c]
    genome_width = 1 - pp["leftmargin"] - pp["rightmargin"] - pp["ideogramlateralmargin"] * n
    mid_y = _ideogram_mid_single_line(pp)

    def ideogram_mid(chrom):
        chrom = np.atleast_1d(np.asarray(chrom, dtype=object))
        return np.full(len(chrom), mid_y, dtype=float)

    def coord_change(chrom=None, x=None, y=None, data_panel=1):
        if plot_type == 4 and (data_panel is None or data_panel == 2):
            data_panel = 1
        if plot_type == 5 and (data_panel is None or data_panel == 1):
            data_panel = 2
        if plot_type == 7 and data_panel in (None, 1, 2):
            data_panel = "ideogram"

        x_plot = None
        y_plot = None

        if x is not None and chrom is not None:
            x = np.asarray(x, dtype=float)
            chrom_arr = np.asarray(chrom, dtype=object)
            idx = np.array([chr_idx[c] for c in chrom_arr], dtype=float)
            starts = np.array([chr_starts[c] for c in chrom_arr], dtype=float)
            prevs = np.array([prev[c] for c in chrom_arr], dtype=float)
            x_plot = (pp["leftmargin"]
                      + pp["ideogramlateralmargin"] * idx
                      + ((prevs + x - starts) / all_chr_len) * genome_width)

        if chrom is not None:
            chrom_arr = np.asarray(chrom, dtype=object)
            mid = np.full(len(chrom_arr), mid_y, dtype=float)
            if y is None:
                y_plot = mid
            else:
                y_arr = np.asarray(y, dtype=float)
                # plot_type 5 inverts y to keep "standard" orientation in panel 2
                if plot_type == 5:
                    y_arr = pp["data2max"] - (y_arr - pp["data2min"])
                y_plot = _scale_y(y_arr, data_panel, pp, mid)

        return x_plot, y_plot

    ylim = (0, pp["bottommargin"] + chr_height + pp["topmargin"])
    return CoordChange(coord_change, ideogram_mid, chr_height, (0.0, 1.0), ylim)
