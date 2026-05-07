"""Ideogram and chromosome-annotation primitives.

Ports of:
  R/kpAddCytobands.R
  R/kpAddCytobandsAsLine.R
  R/kpAddCytobandLabels.R
  R/kpAddChromosomeNames.R
  R/kpAddChromosomeSeparators.R
  R/kpAddBaseNumbers.R
  R/kpAddMainTitle.R
  R/kpAddLabels.R
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

from .colors import get_cytoband_colors
from .ranges import GRanges


def kp_add_cytobands(karyoplot, *, color_table=None, color_schema="circos",
                     border="black", lwd=1.0, clipping=True):
    """Draw cytobands on each chromosome ideogram."""
    cyto = karyoplot.cytobands
    if cyto is None or len(cyto) == 0:
        # synthesize a single gpos50 band per chromosome
        chrom = np.array(karyoplot.chromosomes, dtype=object)
        starts = np.array([karyoplot.chromosome_starts[c] for c in chrom], dtype=np.int64)
        ends = np.array([karyoplot.chromosome_starts[c] + karyoplot.chromosome_lengths[c]
                         for c in chrom], dtype=np.int64)
        cyto = GRanges(chrom, starts, ends,
                       {"name": chrom, "gie_stain": np.array(["gpos50"] * len(chrom), dtype=object)})

    cyto = cyto.filter_chromosomes(karyoplot.chromosomes)
    if len(cyto) == 0:
        return karyoplot

    color_table = get_cytoband_colors(color_table=color_table, color_schema=color_schema)
    if border is None and "border" in color_table:
        border = color_table["border"]

    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    mids = karyoplot.ideogram_mid(cyto.chrom)
    ybottom = mids - pp["ideogramheight"] / 2
    ytop = mids + pp["ideogramheight"] / 2

    xleft, _ = ccf(chrom=cyto.chrom, x=cyto.start, data_panel="ideogram")
    xright, _ = ccf(chrom=cyto.chrom, x=cyto.end, data_panel="ideogram")

    stains = cyto.mcols.get("gie_stain", np.array(["gpos50"] * len(cyto), dtype=object))
    fills = np.array([color_table.get(str(s), "#C8C8C8") for s in stains], dtype=object)

    # Render acen (centromere) bands as triangles, others as rectangles
    is_acen = np.array([str(s) == "acen" for s in stains], dtype=bool)
    rect_mask = ~is_acen
    if rect_mask.any():
        karyoplot.backend.rect(xleft[rect_mask], ybottom[rect_mask],
                               xright[rect_mask], ytop[rect_mask],
                               fill=fills[rect_mask], border=border, lwd=lwd)
    # Acen: draw as filled triangles pointing toward the centromere midpoint per chromosome
    if is_acen.any():
        # group acen by chromosome so we know which side the triangle points
        for chrom_name in np.unique(cyto.chrom[is_acen]):
            mask = is_acen & (cyto.chrom == chrom_name)
            if not mask.any():
                continue
            xs_l = xleft[mask]; xs_r = xright[mask]
            yb = ybottom[mask][0]; yt = ytop[mask][0]
            mid_x = (xs_l.min() + xs_r.max()) / 2
            for i in range(len(xs_l)):
                xl_i, xr_i = float(xs_l[i]), float(xs_r[i])
                if (xl_i + xr_i) / 2 < mid_x:
                    pts_x = [xl_i, xr_i, xl_i]
                    pts_y = [yb, (yb + yt) / 2, yt]
                else:
                    pts_x = [xr_i, xl_i, xr_i]
                    pts_y = [yb, (yb + yt) / 2, yt]
                karyoplot.backend.polygon(pts_x, pts_y,
                                          fill=fills[mask][i], border=border, lwd=lwd)
    return karyoplot


def kp_add_cytobands_as_line(karyoplot, *, color_table=None, color_schema="circos",
                             lwd=3.0):
    """Render the chromosome as a thin horizontal line (no per-band rectangles)."""
    color_table = get_cytoband_colors(color_table=color_table, color_schema=color_schema)
    ccf = karyoplot.coord_change
    for c in karyoplot.chromosomes:
        s = karyoplot.chromosome_starts[c]
        e = s + karyoplot.chromosome_lengths[c]
        xl, _ = ccf(chrom=[c], x=[s], data_panel="ideogram")
        xr, _ = ccf(chrom=[c], x=[e], data_panel="ideogram")
        mid = karyoplot.ideogram_mid([c])[0]
        karyoplot.backend.segments([xl[0]], [mid], [xr[0]], [mid],
                                   color=color_table.get("gpos50", "#888888"), lwd=lwd)
    return karyoplot


def kp_add_cytoband_labels(karyoplot, *, srt: float = 0.0,
                           cex: float = 0.5, color="black"):
    """Draw cytoband names below each band."""
    cyto = karyoplot.cytobands
    if cyto is None or len(cyto) == 0:
        return karyoplot
    cyto = cyto.filter_chromosomes(karyoplot.chromosomes)
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    mids = karyoplot.ideogram_mid(cyto.chrom)
    xs = (cyto.start + cyto.end) / 2
    xp, _ = ccf(chrom=cyto.chrom, x=xs, data_panel="ideogram")
    yp = mids + pp["ideogramheight"] / 2 + 4
    names = cyto.mcols.get("name", np.array([""] * len(cyto), dtype=object))
    karyoplot.backend.text(xp, yp, list(names), color=color, size=10 * cex,
                           halign="center", valign="bottom", angle=srt)
    return karyoplot


def kp_add_chromosome_names(karyoplot, *, chr_names=None, color="black",
                            cex: float = 1.0, srt: float = 0.0):
    """Draw chromosome names to the left of (or above, single-line) each ideogram."""
    chroms = list(karyoplot.chromosomes)
    if chr_names is None:
        labels = chroms
    else:
        labels = list(chr_names)
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    if karyoplot.plot_type in (3, 4, 5, 7):
        # single-line: place names above each chromosome at its midpoint
        for c, lab in zip(chroms, labels):
            s = karyoplot.chromosome_starts[c]
            e = s + karyoplot.chromosome_lengths[c]
            xm, _ = ccf(chrom=[c], x=[(s + e) / 2], data_panel="ideogram")
            mid = karyoplot.ideogram_mid([c])[0]
            yp = mid + pp["ideogramheight"] / 2 + 6
            karyoplot.backend.text([xm[0]], [yp], [lab], color=color,
                                   size=12 * cex, halign="center", valign="bottom",
                                   angle=srt)
    else:
        for c, lab in zip(chroms, labels):
            mid = karyoplot.ideogram_mid([c])[0]
            x = pp["leftmargin"] - 0.005
            karyoplot.backend.text([x], [mid], [lab], color=color,
                                   size=12 * cex, halign="right", valign="middle",
                                   angle=srt)
    return karyoplot


def kp_add_chromosome_separators(karyoplot, *, color="#888888", lwd=1.0,
                                 lty="solid"):
    """Vertical separator between adjacent chromosomes (single-line layouts only)."""
    if karyoplot.plot_type not in (3, 4, 5, 7):
        return karyoplot
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    chroms = karyoplot.chromosomes
    ylim = karyoplot.coord.plot_ylim
    for c0, c1 in zip(chroms[:-1], chroms[1:]):
        e0 = karyoplot.chromosome_starts[c0] + karyoplot.chromosome_lengths[c0]
        s1 = karyoplot.chromosome_starts[c1]
        x0_p, _ = ccf(chrom=[c0], x=[e0], data_panel="ideogram")
        x1_p, _ = ccf(chrom=[c1], x=[s1], data_panel="ideogram")
        x = (float(x0_p[0]) + float(x1_p[0])) / 2
        karyoplot.backend.segments([x], [ylim[0] + pp["bottommargin"]],
                                   [x], [ylim[1] - pp["topmargin"]],
                                   color=color, lwd=lwd, lty=lty)
    return karyoplot


def kp_add_base_numbers(karyoplot, *, tick_dist: int = 20_000_000,
                        tick_len: float = 0.005, add_text: bool = True,
                        label_offset: float = 0.003,
                        color="black", cex: float = 0.6):
    """Tick marks below each ideogram at multiples of ``tick_dist``.

    ``tick_len`` and ``label_offset`` are in figure-fraction units (the same
    coordinate system the coord_change functions return).
    """
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    for c in karyoplot.chromosomes:
        start = karyoplot.chromosome_starts[c]
        end = start + karyoplot.chromosome_lengths[c]
        rem = start % tick_dist
        first = start if rem == 0 else start + (tick_dist - rem)
        ticks = np.arange(first, end + 1, tick_dist, dtype=np.int64)
        if len(ticks) == 0:
            continue
        chr_arr = np.array([c] * len(ticks), dtype=object)
        xp, _ = ccf(chrom=chr_arr, x=ticks, data_panel="ideogram")
        mid = karyoplot.ideogram_mid([c])[0]
        ybot = mid - pp["ideogramheight"] / 2
        karyoplot.backend.segments(xp, np.full(len(xp), ybot),
                                   xp, np.full(len(xp), ybot - tick_len),
                                   color=color, lwd=0.8)
        if add_text:
            labels = [f"{t // 1_000_000} Mb" for t in ticks]
            karyoplot.backend.text(xp,
                                   np.full(len(xp), ybot - tick_len - label_offset),
                                   labels, color=color, size=10 * cex,
                                   halign="center", valign="top")
    return karyoplot


def kp_add_main_title(karyoplot, main: str, *, color="black", cex: float = 1.5):
    pp = karyoplot.plot_params
    ylim = karyoplot.coord.plot_ylim
    yp = ylim[1] - pp["topmargin"] / 2
    karyoplot.backend.text([0.5], [yp], [main], color=color,
                           size=18 * cex, halign="center", valign="middle")
    return karyoplot


def kp_add_labels(karyoplot, label: str, *, data_panel=1, r0=None, r1=None,
                  color="black", cex: float = 1.0, side: str = "left",
                  srt: float = 0.0):
    """Add a label to the left (or right) margin of a data panel.

    In stacked layouts the label is drawn next to every chromosome's panel.
    """
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    from .prepare_params import preprocess_r0_r1
    r0, r1 = preprocess_r0_r1(r0, r1, default=(0.0, 1.0))
    mid_y_data = (r0 + r1) / 2
    if side == "left":
        x = pp["leftmargin"] - 0.02
        halign = "right"
    else:
        x = 1 - pp["rightmargin"] + 0.02
        halign = "left"
    chroms = karyoplot.chromosomes
    _, yp = ccf(chrom=chroms, y=np.full(len(chroms), mid_y_data),
                data_panel=data_panel)
    # Single-line layouts collapse all chromosomes onto one band — dedupe.
    seen = set()
    for y_val in yp:
        key = round(float(y_val), 6)
        if key in seen:
            continue
        seen.add(key)
        karyoplot.backend.text([x], [float(y_val)], [label], color=color,
                               size=12 * cex, halign=halign, valign="middle",
                               angle=srt)
    return karyoplot
