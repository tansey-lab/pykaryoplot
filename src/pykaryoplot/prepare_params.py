"""Parameter normalization helpers (port of prepareParameters2/4 + utils)."""
from __future__ import annotations

from typing import Any, Optional

import numpy as np

from .ranges import GRanges, to_granges


def _panel_range(pp: dict, data_panel) -> tuple[float, float]:
    if data_panel == 1:
        return pp["data1min"], pp["data1max"]
    if data_panel == 2:
        return pp["data2min"], pp["data2max"]
    if data_panel == "ideogram":
        return pp["dataideogrammin"], pp["dataideogrammax"]
    return pp["dataallmin"], pp["dataallmax"]


def preprocess_r0_r1(r0, r1, default=(0.0, 1.0)):
    """Mirror of karyoploteR's r0/r1 handling.

    If ``r1`` is None and ``r0`` is a (r0, r1) sequence or a dict-like with
    ``r0``/``r1`` keys, expand it. Falls back to ``default``.
    """
    if r1 is None:
        if isinstance(r0, dict) and "r0" in r0 and "r1" in r0:
            return float(r0["r0"]), float(r0["r1"])
        if isinstance(r0, (list, tuple, np.ndarray)) and len(r0) == 2:
            return float(r0[0]), float(r0[1])
    if r0 is None:
        r0 = default[0]
    if r1 is None:
        r1 = default[1]
    return float(r0), float(r1)


def _resolve_chr_x_y(data, chrom, x, y):
    """Resolve (chrom, x, y) from a data GRanges or explicit args."""
    if data is not None:
        gr = to_granges(data)
        if chrom is None:
            chrom = gr.chrom
        if x is None:
            x = gr.midpoints
        if y is None:
            for k in ("y", "value"):
                if k in gr.mcols:
                    y = gr.mcols[k]
                    break
    return chrom, x, y


def prepare_parameters_2(karyoplot, *, data=None, chrom=None, x=None, y=None,
                         ymin=None, ymax=None, r0=None, r1=None, data_panel=1):
    """Normalize chr/x/y for point/line-style plots.

    Returns a namespace with ``chrom``, ``x``, ``y_scaled`` (already mapped
    into the data-panel y range), ``filter`` (a boolean mask of valid points),
    and ``original_length``.
    """
    chrom, x, y = _resolve_chr_x_y(data, chrom, x, y)
    if chrom is None:
        raise ValueError("kp_* functions need either `data` or `chrom`")

    chrom = np.asarray(chrom, dtype=object)
    n_orig = len(chrom)

    if x is not None:
        x = np.broadcast_to(np.asarray(x, dtype=float), (n_orig,)).copy()
    if y is not None:
        y = np.broadcast_to(np.asarray(y, dtype=float), (n_orig,)).copy()

    pp = karyoplot.plot_params
    panel_default_min, panel_default_max = _panel_range(pp, data_panel)
    d_min = panel_default_min if ymin is None else ymin
    d_max = panel_default_max if ymax is None else ymax

    r0, r1 = preprocess_r0_r1(r0, r1, default=(panel_default_min, panel_default_max))

    if y is not None and (d_max - d_min) != 0:
        # Map y from [ymin, ymax] -> [0, 1] -> [r0, r1] in panel-canonical coords
        t = (y - d_min) / (d_max - d_min)
        y_scaled = r0 + t * (r1 - r0)
    else:
        y_scaled = y

    # filter to chromosomes that are actually visible
    visible = set(karyoplot.chromosomes)
    chr_filter = np.array([c in visible for c in chrom], dtype=bool)

    return _Params(chrom=chrom, x=x, y=y_scaled,
                   filter=chr_filter, original_length=n_orig,
                   ymin=d_min, ymax=d_max, r0=r0, r1=r1)


def prepare_parameters_4(karyoplot, *, data=None, chrom=None, x0=None, x1=None,
                         y0=None, y1=None, ymin=None, ymax=None,
                         r0=None, r1=None, data_panel=1):
    """Normalize chr/x0/x1/y0/y1 for rectangle-style plots."""
    if data is not None:
        gr = to_granges(data)
        if chrom is None:
            chrom = gr.chrom
        if x0 is None:
            x0 = gr.start
        if x1 is None:
            x1 = gr.end
        if y0 is None and "y0" in gr.mcols:
            y0 = gr.mcols["y0"]
        if y1 is None and "y1" in gr.mcols:
            y1 = gr.mcols["y1"]
    if chrom is None:
        raise ValueError("kp_* functions need either `data` or `chrom`")
    chrom = np.asarray(chrom, dtype=object)
    n_orig = len(chrom)
    if x0 is not None:
        x0 = np.broadcast_to(np.asarray(x0, dtype=float), (n_orig,)).copy()
    if x1 is not None:
        x1 = np.broadcast_to(np.asarray(x1, dtype=float), (n_orig,)).copy()

    pp = karyoplot.plot_params
    panel_default_min, panel_default_max = _panel_range(pp, data_panel)
    d_min = panel_default_min if ymin is None else ymin
    d_max = panel_default_max if ymax is None else ymax

    r0, r1 = preprocess_r0_r1(r0, r1, default=(panel_default_min, panel_default_max))

    def _scale(v):
        if v is None:
            return None
        v = np.broadcast_to(np.asarray(v, dtype=float), (n_orig,)).copy()
        if (d_max - d_min) == 0:
            return v
        t = (v - d_min) / (d_max - d_min)
        return r0 + t * (r1 - r0)

    if y0 is None:
        y0 = np.full(n_orig, r0, dtype=float)
    else:
        y0 = _scale(y0)
    if y1 is None:
        y1 = np.full(n_orig, r1, dtype=float)
    else:
        y1 = _scale(y1)

    visible = set(karyoplot.chromosomes)
    chr_filter = np.array([c in visible for c in chrom], dtype=bool)

    return _Params4(chrom=chrom, x0=x0, x1=x1, y0=y0, y1=y1,
                    filter=chr_filter, original_length=n_orig,
                    ymin=d_min, ymax=d_max, r0=r0, r1=r1)


class _Params:
    __slots__ = ("chrom", "x", "y", "filter", "original_length",
                 "ymin", "ymax", "r0", "r1")

    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


class _Params4(_Params):
    __slots__ = ("chrom", "x0", "x1", "y0", "y1", "filter", "original_length",
                 "ymin", "ymax", "r0", "r1")


def autotrack(current: int, total: int, margin: float = 0.05,
              r0: float = 0.0, r1: float = 1.0) -> tuple[float, float]:
    """Compute (r0, r1) for the i-th of N stacked tracks within [r0,r1]."""
    if current < 1 or current > total:
        raise ValueError("current must be in 1..total")
    span = (r1 - r0) / total
    lo = r0 + (current - 1) * span + margin * span / 2
    hi = r0 + current * span - margin * span / 2
    return lo, hi


def filter_array(arr, mask: np.ndarray, original_length: int):
    """If ``arr`` is a per-element array, subset it; otherwise pass through."""
    if arr is None:
        return None
    if isinstance(arr, (list, tuple, np.ndarray)) and len(arr) == original_length:
        return np.asarray(arr)[mask]
    return arr


def apply_clip(karyoplot, data_panel) -> None:
    """If zoom is active, clip subsequent drawing to the data panel area."""
    if not karyoplot.zoom:
        return
    pp = karyoplot.plot_params
    region = karyoplot.plot_region
    chrom = str(region.chrom[0])
    ccf = karyoplot.coord_change
    xl, _ = ccf(chrom=[chrom], x=[int(region.start[0])], data_panel=data_panel)
    xr, _ = ccf(chrom=[chrom], x=[int(region.end[0])], data_panel=data_panel)
    mid = karyoplot.ideogram_mid([chrom])[0]
    if data_panel == 1:
        yb = mid + pp["ideogramheight"] / 2 + pp["data1inmargin"]
        yt = yb + pp["data1height"]
    elif data_panel == 2:
        yt = mid - pp["ideogramheight"] / 2 - pp["data2inmargin"]
        yb = yt - pp["data2height"]
    elif data_panel == "ideogram":
        yb = mid - pp["ideogramheight"] / 2
        yt = mid + pp["ideogramheight"] / 2
    else:
        return
    karyoplot.backend.push_clip(float(xl[0]), float(yb), float(xr[0]), float(yt))
