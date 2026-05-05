"""Data plotting primitives — kp_points / kp_lines / kp_rect / ..."""
from __future__ import annotations

from typing import Optional

import numpy as np

from .colors import col_by_value
from .prepare_params import (apply_clip, prepare_parameters_2,
                             prepare_parameters_4, preprocess_r0_r1, filter_array)


# --- Helpers --------------------------------------------------------------

def _maybe_clip(karyoplot, clipping, data_panel):
    if clipping and karyoplot.zoom:
        apply_clip(karyoplot, data_panel)
        return True
    return False


def _maybe_unclip(karyoplot, did_clip):
    if did_clip:
        karyoplot.backend.pop_clip()


# --- Point/line family ----------------------------------------------------

def kp_points(karyoplot, *, data=None, chrom=None, x=None, y=None,
              ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
              clipping=True, color="black", size: float = 2.5, pch: int = 16):
    p = prepare_parameters_2(karyoplot, data=data, chrom=chrom, x=x, y=y,
                             ymin=ymin, ymax=ymax, r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    xp, yp = ccf(chrom=p.chrom, x=p.x, y=p.y, data_panel=data_panel)
    color_arr = filter_array(color, p.filter, p.original_length)
    did = _maybe_clip(karyoplot, clipping, data_panel)
    karyoplot.backend.points(xp[p.filter], yp[p.filter],
                             color=color_arr if isinstance(color_arr, np.ndarray) else color,
                             size=size, pch=pch)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_lines(karyoplot, *, data=None, chrom=None, x=None, y=None,
             ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
             clipping=True, color="black", lwd: float = 1.0, lty="solid"):
    p = prepare_parameters_2(karyoplot, data=data, chrom=chrom, x=x, y=y,
                             ymin=ymin, ymax=ymax, r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    did = _maybe_clip(karyoplot, clipping, data_panel)
    # Plot one line per chromosome
    for c in np.unique(p.chrom[p.filter]):
        m = (p.chrom == c) & p.filter
        if m.sum() < 2:
            continue
        xs = p.x[m]
        ys = p.y[m]
        order = np.argsort(xs)
        xp, yp = ccf(chrom=p.chrom[m][order], x=xs[order], y=ys[order],
                     data_panel=data_panel)
        karyoplot.backend.lines(xp, yp, color=color, lwd=lwd, lty=lty)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_segments(karyoplot, *, data=None, chrom=None, x0=None, x1=None,
                y0=None, y1=None, ymin=None, ymax=None, data_panel=1,
                r0=None, r1=None, clipping=True, color="black",
                lwd: float = 1.0, lty="solid"):
    p = prepare_parameters_4(karyoplot, data=data, chrom=chrom, x0=x0, x1=x1,
                             y0=y0, y1=y1, ymin=ymin, ymax=ymax,
                             r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    x0p, y0p = ccf(chrom=p.chrom, x=p.x0, y=p.y0, data_panel=data_panel)
    x1p, y1p = ccf(chrom=p.chrom, x=p.x1, y=p.y1, data_panel=data_panel)
    did = _maybe_clip(karyoplot, clipping, data_panel)
    karyoplot.backend.segments(x0p[p.filter], y0p[p.filter],
                               x1p[p.filter], y1p[p.filter],
                               color=color, lwd=lwd, lty=lty)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_arrows(karyoplot, *, data=None, chrom=None, x0=None, x1=None,
              y0=None, y1=None, ymin=None, ymax=None, data_panel=1,
              r0=None, r1=None, clipping=True, color="black",
              lwd: float = 1.0, head_length: float = 6.0):
    p = prepare_parameters_4(karyoplot, data=data, chrom=chrom, x0=x0, x1=x1,
                             y0=y0, y1=y1, ymin=ymin, ymax=ymax,
                             r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    x0p, y0p = ccf(chrom=p.chrom, x=p.x0, y=p.y0, data_panel=data_panel)
    x1p, y1p = ccf(chrom=p.chrom, x=p.x1, y=p.y1, data_panel=data_panel)
    did = _maybe_clip(karyoplot, clipping, data_panel)
    karyoplot.backend.arrows(x0p[p.filter], y0p[p.filter],
                             x1p[p.filter], y1p[p.filter],
                             color=color, lwd=lwd, head_length=head_length)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_text(karyoplot, *, labels, data=None, chrom=None, x=None, y=None,
            ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
            clipping=True, color="black", cex: float = 1.0, srt: float = 0.0,
            halign="center", valign="middle"):
    p = prepare_parameters_2(karyoplot, data=data, chrom=chrom, x=x, y=y,
                             ymin=ymin, ymax=ymax, r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    xp, yp = ccf(chrom=p.chrom, x=p.x, y=p.y, data_panel=data_panel)
    labs = list(labels) if not isinstance(labels, str) else [labels] * p.original_length
    labs = np.asarray(labs)[p.filter]
    did = _maybe_clip(karyoplot, clipping, data_panel)
    karyoplot.backend.text(xp[p.filter], yp[p.filter], list(labs),
                           color=color, size=12 * cex, halign=halign,
                           valign=valign, angle=srt)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_abline(karyoplot, *, h=None, v=None, chrom=None, ymin=None, ymax=None,
              data_panel=1, r0=None, r1=None, color="black",
              lwd: float = 1.0, lty="solid"):
    """Horizontal/vertical reference lines.

    ``h``: data-panel y values (broadcast across all visible chromosomes
    unless ``chrom`` is given).
    ``v``: genomic positions in bp, requires ``chrom``.
    """
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    if h is not None:
        h = np.atleast_1d(np.asarray(h, dtype=float))
        chroms = karyoplot.chromosomes if chrom is None else list(np.atleast_1d(chrom))
        for c in chroms:
            for hv in h:
                _, yp = ccf(chrom=[c], y=[hv], data_panel=data_panel)
                start = karyoplot.chromosome_starts[c]
                end = start + karyoplot.chromosome_lengths[c]
                xl, _ = ccf(chrom=[c], x=[start], data_panel=data_panel)
                xr, _ = ccf(chrom=[c], x=[end], data_panel=data_panel)
                karyoplot.backend.segments([xl[0]], [yp[0]], [xr[0]], [yp[0]],
                                           color=color, lwd=lwd, lty=lty)
    if v is not None:
        if chrom is None:
            raise ValueError("kp_abline(v=...) requires chrom=")
        v = np.atleast_1d(np.asarray(v, dtype=float))
        chrom_arr = np.atleast_1d(np.asarray(chrom, dtype=object))
        if len(chrom_arr) == 1 and len(v) > 1:
            chrom_arr = np.repeat(chrom_arr, len(v))
        for c, vv in zip(chrom_arr, v):
            xp, _ = ccf(chrom=[c], x=[vv], data_panel=data_panel)
            mid = karyoplot.ideogram_mid([c])[0]
            yt = mid + pp["ideogramheight"] / 2 + pp["data1inmargin"] + pp["data1height"]
            yb = mid - pp["ideogramheight"] / 2 - pp["data2inmargin"] - pp["data2height"]
            karyoplot.backend.segments([xp[0]], [yb], [xp[0]], [yt],
                                       color=color, lwd=lwd, lty=lty)
    return karyoplot


# --- Rectangle family -----------------------------------------------------

def kp_rect(karyoplot, *, data=None, chrom=None, x0=None, x1=None,
            y0=None, y1=None, ymin=None, ymax=None, data_panel=1,
            r0=None, r1=None, clipping=True, color="#AAAAAA", border="black",
            lwd: float = 1.0):
    p = prepare_parameters_4(karyoplot, data=data, chrom=chrom, x0=x0, x1=x1,
                             y0=y0, y1=y1, ymin=ymin, ymax=ymax,
                             r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    xl, yb = ccf(chrom=p.chrom, x=p.x0, y=p.y0, data_panel=data_panel)
    xr, yt = ccf(chrom=p.chrom, x=p.x1, y=p.y1, data_panel=data_panel)
    did = _maybe_clip(karyoplot, clipping, data_panel)
    karyoplot.backend.rect(xl[p.filter], yb[p.filter], xr[p.filter], yt[p.filter],
                           fill=color, border=border, lwd=lwd)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_bars(karyoplot, *, data=None, chrom=None, x0=None, x1=None, y0=None, y1=None,
            ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
            clipping=True, color="#7777CC", border=None):
    return kp_rect(karyoplot, data=data, chrom=chrom, x0=x0, x1=x1,
                   y0=y0, y1=y1, ymin=ymin, ymax=ymax, data_panel=data_panel,
                   r0=r0, r1=r1, clipping=clipping, color=color, border=border)


def kp_polygon(karyoplot, *, data=None, chrom=None, x=None, y=None,
               ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
               clipping=True, color="#AAAAAA", border="black", lwd: float = 1.0):
    p = prepare_parameters_2(karyoplot, data=data, chrom=chrom, x=x, y=y,
                             ymin=ymin, ymax=ymax, r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    did = _maybe_clip(karyoplot, clipping, data_panel)
    for c in np.unique(p.chrom[p.filter]):
        m = (p.chrom == c) & p.filter
        xp, yp = ccf(chrom=p.chrom[m], x=p.x[m], y=p.y[m], data_panel=data_panel)
        karyoplot.backend.polygon(xp, yp, fill=color, border=border, lwd=lwd)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_area(karyoplot, *, data=None, chrom=None, x=None, y=None,
            base_y: float = 0.0, ymin=None, ymax=None, data_panel=1,
            r0=None, r1=None, clipping=True, color="#88CC8888", border="#226622",
            lwd: float = 1.0):
    """Filled area under a line."""
    p = prepare_parameters_2(karyoplot, data=data, chrom=chrom, x=x, y=y,
                             ymin=ymin, ymax=ymax, r0=r0, r1=r1, data_panel=data_panel)
    ccf = karyoplot.coord_change
    did = _maybe_clip(karyoplot, clipping, data_panel)
    for c in np.unique(p.chrom[p.filter]):
        m = (p.chrom == c) & p.filter
        if m.sum() < 2:
            continue
        order = np.argsort(p.x[m])
        xs = p.x[m][order]; ys = p.y[m][order]
        chrs = p.chrom[m][order]
        xp, yp = ccf(chrom=chrs, x=xs, y=ys, data_panel=data_panel)
        # close to base_y
        _, by = ccf(chrom=[c, c], y=[base_y, base_y], data_panel=data_panel)
        poly_x = np.concatenate([xp, xp[::-1]])
        poly_y = np.concatenate([yp, np.full(len(yp), float(by[0]))])
        karyoplot.backend.polygon(poly_x, poly_y, fill=color, border=border, lwd=lwd)
    _maybe_unclip(karyoplot, did)
    return karyoplot


def kp_heatmap(karyoplot, *, data=None, chrom=None, x0=None, x1=None, y=None,
               colors=("#FFFFFF", "#000000"), ymin=None, ymax=None,
               data_panel=1, r0=None, r1=None, clipping=True):
    """Heatmap rectangles colored by ``y`` value."""
    if y is None and data is not None:
        from .ranges import to_granges
        gr = to_granges(data)
        for k in ("y", "value"):
            if k in gr.mcols:
                y = gr.mcols[k]; break
    p = prepare_parameters_4(karyoplot, data=data, chrom=chrom, x0=x0, x1=x1,
                             y0=None, y1=None, ymin=ymin, ymax=ymax,
                             r0=r0, r1=r1, data_panel=data_panel)
    if y is None:
        raise ValueError("kp_heatmap needs y values")
    y = np.asarray(y, dtype=float)
    fills = col_by_value(y, colors=colors,
                         vmin=ymin if ymin is not None else float(np.nanmin(y)),
                         vmax=ymax if ymax is not None else float(np.nanmax(y)))
    ccf = karyoplot.coord_change
    xl, yb = ccf(chrom=p.chrom, x=p.x0, y=p.y0, data_panel=data_panel)
    xr, yt = ccf(chrom=p.chrom, x=p.x1, y=p.y1, data_panel=data_panel)
    did = _maybe_clip(karyoplot, clipping, data_panel)
    karyoplot.backend.rect(xl[p.filter], yb[p.filter], xr[p.filter], yt[p.filter],
                           fill=fills[p.filter], border=None)
    _maybe_unclip(karyoplot, did)
    return karyoplot


# --- Background / axis / labels -------------------------------------------

def kp_data_background(karyoplot, *, color="#EEEEEE", data_panel=1,
                       r0=None, r1=None):
    pp = karyoplot.plot_params
    ccf = karyoplot.coord_change
    r0, r1 = preprocess_r0_r1(r0, r1, default=(0.0, 1.0))
    if data_panel == 1:
        d_min, d_max = pp["data1min"], pp["data1max"]
    else:
        d_min, d_max = pp["data2min"], pp["data2max"]
    y_lo = d_min + (r0 - 0.0) * (d_max - d_min)
    y_hi = d_min + (r1 - 0.0) * (d_max - d_min)
    for c in karyoplot.chromosomes:
        s = karyoplot.chromosome_starts[c]
        e = s + karyoplot.chromosome_lengths[c]
        xl, ylo = ccf(chrom=[c], x=[s], y=[y_lo], data_panel=data_panel)
        xr, yhi = ccf(chrom=[c], x=[e], y=[y_hi], data_panel=data_panel)
        karyoplot.backend.rect([xl[0]], [ylo[0]], [xr[0]], [yhi[0]],
                               fill=color, border=None)
    return karyoplot


def kp_axis(karyoplot, *, side: int = 1, ymin=None, ymax=None,
            tick_pos=None, n_ticks: int = 5, data_panel=1,
            r0=None, r1=None, color="black", cex: float = 0.6):
    """Vertical axis on the left (side=1) or right (side=2) of each ideogram."""
    pp = karyoplot.plot_params
    if data_panel == 1:
        d_min = pp["data1min"] if ymin is None else ymin
        d_max = pp["data1max"] if ymax is None else ymax
    else:
        d_min = pp["data2min"] if ymin is None else ymin
        d_max = pp["data2max"] if ymax is None else ymax
    if tick_pos is None:
        tick_pos = np.linspace(d_min, d_max, n_ticks)
    tick_pos = np.asarray(tick_pos, dtype=float)
    ccf = karyoplot.coord_change
    for c in karyoplot.chromosomes:
        s = karyoplot.chromosome_starts[c]
        e = s + karyoplot.chromosome_lengths[c]
        x_anchor = s if side == 1 else e
        ha = "right" if side == 1 else "left"
        for tp in tick_pos:
            xp, yp = ccf(chrom=[c], x=[x_anchor], y=[tp], data_panel=data_panel)
            karyoplot.backend.text([float(xp[0]) + (-0.005 if side == 1 else 0.005)],
                                   [float(yp[0])], [f"{tp:g}"],
                                   color=color, size=10 * cex,
                                   halign=ha, valign="middle")
    return karyoplot


def kp_add_color_rect(karyoplot, *, color="black", data_panel=1, r0=None, r1=None,
                      side: str = "left", width: float = 0.005):
    """Color stripe to the side of a data panel (track marker)."""
    pp = karyoplot.plot_params
    r0, r1 = preprocess_r0_r1(r0, r1, default=(0.0, 1.0))
    ccf = karyoplot.coord_change
    if data_panel == 1:
        d_min, d_max = pp["data1min"], pp["data1max"]
    else:
        d_min, d_max = pp["data2min"], pp["data2max"]
    y_lo = d_min + r0 * (d_max - d_min)
    y_hi = d_min + r1 * (d_max - d_min)
    chrom_ref = karyoplot.chromosomes[0]
    _, ylo = ccf(chrom=[chrom_ref], y=[y_lo], data_panel=data_panel)
    _, yhi = ccf(chrom=[chrom_ref], y=[y_hi], data_panel=data_panel)
    if side == "left":
        x0 = pp["leftmargin"] - width
        x1 = pp["leftmargin"]
    else:
        x0 = 1 - pp["rightmargin"]
        x1 = 1 - pp["rightmargin"] + width
    karyoplot.backend.rect([x0], [float(ylo[0])], [x1], [float(yhi[0])],
                           fill=color, border=None)
    return karyoplot
