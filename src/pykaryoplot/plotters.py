"""High-level plotters with no heavy bioinformatics dependencies."""
from __future__ import annotations

import numpy as np

from .primitives import (kp_bars, kp_lines, kp_points, kp_polygon, kp_rect,
                         kp_segments)
from .ranges import GRanges, to_granges


def kp_plot_regions(karyoplot, regions, *, data_panel=1, r0=None, r1=None,
                    color="#AAAAFF", border="#3333AA", num_layers: int = 1,
                    layer_margin: float = 0.05, avoid_overlapping: bool = True):
    """Draw genomic regions as stacked rectangles (single-layer by default)."""
    gr = to_granges(regions)
    if num_layers == 1 and not avoid_overlapping:
        return kp_rect(karyoplot, data=gr, y0=0.1, y1=0.9, ymin=0, ymax=1,
                       data_panel=data_panel, r0=r0, r1=r1, color=color, border=border)
    # Greedy layering: per chromosome, sort by start; place into first track that ends before this start
    chr_layers = {}
    layer_assignment = np.zeros(len(gr), dtype=int)
    for c in np.unique(gr.chrom):
        idx = np.where(gr.chrom == c)[0]
        order = np.argsort(gr.start[idx])
        idx = idx[order]
        layer_ends = []
        for i in idx:
            placed = False
            for li, end in enumerate(layer_ends):
                if gr.start[i] > end:
                    layer_assignment[i] = li
                    layer_ends[li] = gr.end[i]
                    placed = True
                    break
            if not placed:
                layer_assignment[i] = len(layer_ends)
                layer_ends.append(gr.end[i])
        chr_layers[c] = len(layer_ends)
    max_layers = max(chr_layers.values()) if chr_layers else 1
    if num_layers < max_layers:
        num_layers = max_layers
    span = 1.0 / num_layers
    margin = layer_margin * span
    y0 = np.array([la * span + margin for la in layer_assignment])
    y1 = np.array([(la + 1) * span - margin for la in layer_assignment])
    return kp_rect(karyoplot, data=gr, y0=y0, y1=y1, ymin=0, ymax=1,
                   data_panel=data_panel, r0=r0, r1=r1, color=color, border=border)


def kp_plot_density(karyoplot, regions, *, window_size: int = 1_000_000,
                    data_panel=1, r0=None, r1=None, color="#3366AA",
                    border="#112244"):
    """Density of regions in fixed-width windows along each chromosome."""
    gr = to_granges(regions)
    bars_chr, bars_x0, bars_x1, bars_y = [], [], [], []
    max_count = 0
    for c in karyoplot.chromosomes:
        s0 = karyoplot.chromosome_starts[c]
        e0 = s0 + karyoplot.chromosome_lengths[c]
        starts = np.arange(s0, e0, window_size, dtype=np.int64)
        ends = np.minimum(starts + window_size, e0)
        # count regions whose midpoint falls in window
        m = gr.chrom == c
        if not m.any():
            counts = np.zeros(len(starts), dtype=np.int64)
        else:
            mids = (gr.start[m] + gr.end[m]) / 2
            counts = np.histogram(mids, bins=np.append(starts, ends[-1]))[0]
        if counts.max(initial=0) > max_count:
            max_count = int(counts.max(initial=0))
        bars_chr.extend([c] * len(starts))
        bars_x0.extend(starts.tolist())
        bars_x1.extend(ends.tolist())
        bars_y.extend(counts.tolist())
    if not bars_chr:
        return karyoplot
    y_arr = np.array(bars_y, dtype=float)
    return kp_bars(karyoplot, chrom=np.array(bars_chr, dtype=object),
                   x0=np.array(bars_x0), x1=np.array(bars_x1),
                   y0=np.zeros_like(y_arr), y1=y_arr,
                   ymin=0, ymax=max(1, max_count),
                   data_panel=data_panel, r0=r0, r1=r1, color=color, border=border)


def kp_plot_coverage(karyoplot, regions, *, data_panel=1, r0=None, r1=None,
                     color="#5577CC"):
    """Per-base coverage of overlapping regions, drawn as a polygon area."""
    gr = to_granges(regions)
    max_cov = 1
    chunks = []  # (chrom, xs, ys)
    for c in karyoplot.chromosomes:
        m = gr.chrom == c
        if not m.any():
            continue
        starts = gr.start[m]
        ends = gr.end[m]
        # event-list sweep
        events = np.concatenate([np.column_stack([starts, np.ones_like(starts)]),
                                 np.column_stack([ends + 1, -np.ones_like(ends)])])
        events = events[np.argsort(events[:, 0])]
        xs, ys = [], []
        cov = 0
        last_x = None
        for x, d in events:
            if last_x is not None and x != last_x:
                xs.append(last_x); ys.append(cov)
                xs.append(x);      ys.append(cov)
            cov += int(d)
            last_x = x
        if cov != 0 and last_x is not None:
            xs.append(last_x); ys.append(cov)
        if ys:
            max_cov = max(max_cov, max(ys))
        chunks.append((c, np.array(xs, dtype=float), np.array(ys, dtype=float)))

    pp_data1 = (karyoplot.plot_params["data1min"], karyoplot.plot_params["data1max"]) \
        if data_panel == 1 else (karyoplot.plot_params["data2min"], karyoplot.plot_params["data2max"])
    from .prepare_params import preprocess_r0_r1
    r0v, r1v = preprocess_r0_r1(r0, r1, default=(0.0, 1.0))
    ccf = karyoplot.coord_change
    for c, xs, ys in chunks:
        if len(xs) == 0:
            continue
        # scale ys 0..max_cov into [r0v,r1v] within data panel
        y_scaled = pp_data1[0] + (ys / max_cov) * (r1v - r0v) * (pp_data1[1] - pp_data1[0]) + r0v * (pp_data1[1] - pp_data1[0])
        chrs = np.array([c] * len(xs), dtype=object)
        xp, yp = ccf(chrom=chrs, x=xs, y=y_scaled, data_panel=data_panel)
        # close polygon to baseline
        _, baseline = ccf(chrom=[c], y=[pp_data1[0] + r0v * (pp_data1[1] - pp_data1[0])],
                          data_panel=data_panel)
        poly_x = np.concatenate([[xp[0]], xp, [xp[-1]]])
        poly_y = np.concatenate([[float(baseline[0])], yp, [float(baseline[0])]])
        karyoplot.backend.polygon(poly_x, poly_y, fill=color, border=color, lwd=1)
    return karyoplot


def kp_plot_markers(karyoplot, *, data=None, chrom=None, x=None, labels,
                    y=0.7, marker_color="black", text_color="black",
                    line_color="black", cex: float = 0.8, data_panel=1,
                    r0=None, r1=None):
    """Vertical needles with text labels (a simple, non-overlap-aware version)."""
    if data is not None:
        gr = to_granges(data)
        if chrom is None: chrom = gr.chrom
        if x is None:     x = gr.midpoints
        if labels is None and "name" in gr.mcols:
            labels = gr.mcols["name"]
    chrom = np.asarray(chrom, dtype=object)
    x = np.asarray(x, dtype=float)
    labels = list(labels)
    n = len(chrom)
    y_arr = np.broadcast_to(np.asarray(y, dtype=float), (n,))
    # Stem: from base (y=0) to y
    kp_segments(karyoplot, chrom=chrom, x0=x, x1=x, y0=np.zeros(n), y1=y_arr,
                ymin=0, ymax=1, data_panel=data_panel, r0=r0, r1=r1,
                color=line_color, lwd=1.0)
    from .primitives import kp_points, kp_text
    kp_points(karyoplot, chrom=chrom, x=x, y=y_arr, ymin=0, ymax=1,
              data_panel=data_panel, r0=r0, r1=r1, color=marker_color,
              size=3, pch=16)
    kp_text(karyoplot, labels=labels, chrom=chrom, x=x, y=y_arr + 0.05,
            ymin=0, ymax=1, data_panel=data_panel, r0=r0, r1=r1,
            color=text_color, cex=cex, halign="center", valign="bottom")
    return karyoplot


def kp_plot_2lines(karyoplot, *, data=None, chrom=None, x=None, y1=None, y2=None,
                   ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
                   color1="#3366AA", color2="#AA3333", lwd: float = 1.0):
    """Plot two lines on the same panel."""
    kp_lines(karyoplot, data=data, chrom=chrom, x=x, y=y1,
             ymin=ymin, ymax=ymax, data_panel=data_panel, r0=r0, r1=r1,
             color=color1, lwd=lwd)
    kp_lines(karyoplot, data=data, chrom=chrom, x=x, y=y2,
             ymin=ymin, ymax=ymax, data_panel=data_panel, r0=r0, r1=r1,
             color=color2, lwd=lwd)
    return karyoplot


def kp_plot_rainfall(karyoplot, variants, *, data_panel=1, r0=None, r1=None,
                     color="#222222", size: float = 1.5):
    """Plot log10 distance to previous variant against position."""
    gr = to_granges(variants)
    chrom_out, x_out, y_out = [], [], []
    max_log = 0
    for c in karyoplot.chromosomes:
        m = gr.chrom == c
        if m.sum() < 2:
            continue
        pos = np.sort(gr.midpoints[m])
        d = np.diff(pos)
        d = np.maximum(d, 1)
        logd = np.log10(d)
        max_log = max(max_log, float(logd.max()))
        chrom_out.extend([c] * len(d))
        x_out.extend(pos[1:].tolist())
        y_out.extend(logd.tolist())
    if not chrom_out:
        return karyoplot
    return kp_points(karyoplot, chrom=np.array(chrom_out, dtype=object),
                     x=np.array(x_out), y=np.array(y_out),
                     ymin=0, ymax=max(1, max_log),
                     data_panel=data_panel, r0=r0, r1=r1,
                     color=color, size=size)


def kp_plot_manhattan(karyoplot, *, data=None, chrom=None, x=None, pval=None,
                      data_panel=1, r0=None, r1=None, color=None,
                      threshold: float = 5e-8, threshold_color="red",
                      size: float = 1.5, ymax=None):
    """Manhattan plot: -log10(pval) per locus, alternating colors per chromosome."""
    if data is not None:
        gr = to_granges(data)
        if chrom is None: chrom = gr.chrom
        if x is None:     x = gr.midpoints
        if pval is None:
            for k in ("pval", "p_value", "p", "y"):
                if k in gr.mcols:
                    pval = gr.mcols[k]; break
    chrom = np.asarray(chrom, dtype=object)
    x = np.asarray(x, dtype=float)
    pval = np.asarray(pval, dtype=float)
    y = -np.log10(np.maximum(pval, 1e-300))
    if ymax is None:
        ymax = float(np.nanmax(y)) * 1.05 if len(y) else 8.0
    if color is None:
        from .colors import col_by_chr
        color = col_by_chr(chrom, ["#3366AA", "#AA9933"])
    kp_points(karyoplot, chrom=chrom, x=x, y=y, ymin=0, ymax=ymax,
              data_panel=data_panel, r0=r0, r1=r1, color=color, size=size, pch=16)
    if threshold is not None:
        from .primitives import kp_abline
        kp_abline(karyoplot, h=-np.log10(threshold), ymin=0, ymax=ymax,
                  data_panel=data_panel, r0=r0, r1=r1,
                  color=threshold_color, lwd=1.0, lty="dashed")
    return karyoplot


def kp_plot_ribbon(karyoplot, *, data=None, chrom=None, x=None, y_lo=None, y_hi=None,
                   ymin=None, ymax=None, data_panel=1, r0=None, r1=None,
                   color="#88AAEE88", border="#3355BB"):
    """Plot a filled ribbon between two y-tracks (per chromosome)."""
    if data is not None:
        gr = to_granges(data)
        if chrom is None: chrom = gr.chrom
        if x is None:     x = gr.midpoints
        if y_lo is None and "y_lo" in gr.mcols: y_lo = gr.mcols["y_lo"]
        if y_hi is None and "y_hi" in gr.mcols: y_hi = gr.mcols["y_hi"]
    chrom = np.asarray(chrom, dtype=object)
    x = np.asarray(x, dtype=float)
    y_lo = np.asarray(y_lo, dtype=float)
    y_hi = np.asarray(y_hi, dtype=float)
    ccf = karyoplot.coord_change
    for c in np.unique(chrom):
        m = chrom == c
        order = np.argsort(x[m])
        xs = x[m][order]
        lo = y_lo[m][order]; hi = y_hi[m][order]
        chrs = chrom[m][order]
        xp1, yp_lo = ccf(chrom=chrs, x=xs, y=lo, data_panel=data_panel)
        xp2, yp_hi = ccf(chrom=chrs, x=xs, y=hi, data_panel=data_panel)
        poly_x = np.concatenate([xp1, xp2[::-1]])
        poly_y = np.concatenate([yp_lo, yp_hi[::-1]])
        karyoplot.backend.polygon(poly_x, poly_y, fill=color, border=border, lwd=1.0)
    return karyoplot
