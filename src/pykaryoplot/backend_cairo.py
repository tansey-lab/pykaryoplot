"""Cairo rendering backend.

Plot coordinates use a base-R-like convention:

* origin ``(0, 0)`` is at the **lower-left** of the canvas;
* ``y`` grows upward;
* x-range is ``[0, 1]`` (relative width — so all ``x`` values produced by the
  coordinate-change function are simply percentages of the plot width);
* y-range is the ``ylim`` returned from :func:`build_coord_change` and is
  measured in the karyoploteR plot-param "units" (heights, margins, etc.,
  expressed in those raw numbers).

Cairo's native y-axis grows downward, so the backend flips internally with
``ctx.transform``.
"""
from __future__ import annotations

import io
import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import cairo
import numpy as np

from .colors import to_rgba


@dataclass
class CairoBackend:
    width_px: int
    height_px: int
    surface: cairo.Surface
    ctx: cairo.Context
    xlim: tuple[float, float]
    ylim: tuple[float, float]
    surface_kind: str = "png"  # png | svg | pdf
    _buffer: io.BytesIO | None = None

    @classmethod
    def create(cls, xlim, ylim, width_px=1200, height_px=800, kind="png") -> "CairoBackend":
        kind = kind.lower()
        buf: io.BytesIO | None = None
        if kind == "png":
            surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, width_px, height_px)
        elif kind == "svg":
            buf = io.BytesIO()
            surface = cairo.SVGSurface(buf, width_px, height_px)
        elif kind == "pdf":
            buf = io.BytesIO()
            surface = cairo.PDFSurface(buf, width_px, height_px)
        else:
            raise ValueError(f"Unknown surface kind {kind!r}")
        ctx = cairo.Context(surface)
        # White background
        ctx.set_source_rgb(1, 1, 1)
        ctx.paint()
        be = cls(width_px, height_px, surface, ctx, xlim, ylim, kind, buf)
        return be

    # ------------------------------------------------------------------
    # Coordinate mapping
    # ------------------------------------------------------------------
    def _map_x(self, x):
        x = np.asarray(x, dtype=float)
        return (x - self.xlim[0]) / (self.xlim[1] - self.xlim[0]) * self.width_px

    def _map_y(self, y):
        y = np.asarray(y, dtype=float)
        # invert because cairo y grows down
        return self.height_px - (y - self.ylim[0]) / (self.ylim[1] - self.ylim[0]) * self.height_px

    def _set_color(self, color, alpha: float = 1.0) -> None:
        r, g, b, a = to_rgba(color)
        self.ctx.set_source_rgba(r, g, b, a * alpha)

    # ------------------------------------------------------------------
    # Primitives
    # ------------------------------------------------------------------
    def rect(self, xleft, ybottom, xright, ytop, *, fill=None, border="black",
             lwd: float = 1.0):
        xleft = self._map_x(xleft)
        xright = self._map_x(xright)
        ytop = self._map_y(ytop)
        ybottom = self._map_y(ybottom)
        n = len(np.atleast_1d(xleft))
        for i in range(n):
            xl = float(np.atleast_1d(xleft)[i])
            xr = float(np.atleast_1d(xright)[i])
            yt = float(np.atleast_1d(ytop)[i])
            yb = float(np.atleast_1d(ybottom)[i])
            self.ctx.rectangle(xl, yt, xr - xl, yb - yt)
            if fill is not None:
                col = fill[i] if isinstance(fill, (list, np.ndarray)) else fill
                self._set_color(col)
                self.ctx.fill_preserve()
            if border is not None:
                col = border[i] if isinstance(border, (list, np.ndarray)) else border
                self._set_color(col)
                self.ctx.set_line_width(lwd)
                self.ctx.stroke()
            else:
                self.ctx.new_path()

    def lines(self, x, y, *, color="black", lwd: float = 1.0, lty="solid"):
        xs = self._map_x(x)
        ys = self._map_y(y)
        if len(xs) < 2:
            return
        self._set_color(color)
        self.ctx.set_line_width(lwd)
        self._set_dash(lty)
        self.ctx.move_to(float(xs[0]), float(ys[0]))
        for i in range(1, len(xs)):
            self.ctx.line_to(float(xs[i]), float(ys[i]))
        self.ctx.stroke()
        self.ctx.set_dash([])

    def segments(self, x0, y0, x1, y1, *, color="black", lwd: float = 1.0, lty="solid"):
        x0a = np.atleast_1d(self._map_x(x0))
        x1a = np.atleast_1d(self._map_x(x1))
        y0a = np.atleast_1d(self._map_y(y0))
        y1a = np.atleast_1d(self._map_y(y1))
        n = max(len(x0a), len(x1a), len(y0a), len(y1a))
        self.ctx.set_line_width(lwd)
        self._set_dash(lty)
        for i in range(n):
            col = color[i] if isinstance(color, (list, np.ndarray)) else color
            self._set_color(col)
            self.ctx.move_to(float(x0a[i % len(x0a)]), float(y0a[i % len(y0a)]))
            self.ctx.line_to(float(x1a[i % len(x1a)]), float(y1a[i % len(y1a)]))
            self.ctx.stroke()
        self.ctx.set_dash([])

    def points(self, x, y, *, color="black", size: float = 3.0, pch: int = 16):
        xs = np.atleast_1d(self._map_x(x))
        ys = np.atleast_1d(self._map_y(y))
        for i in range(len(xs)):
            col = color[i] if isinstance(color, (list, np.ndarray)) else color
            self._set_color(col)
            cx, cy = float(xs[i]), float(ys[i])
            r = size
            if pch in (16, 19, 20):  # solid circle
                self.ctx.arc(cx, cy, r, 0, 2 * math.pi)
                self.ctx.fill()
            elif pch == 1:  # open circle
                self.ctx.arc(cx, cy, r, 0, 2 * math.pi)
                self.ctx.set_line_width(1.0)
                self.ctx.stroke()
            elif pch in (15, 22):  # solid square
                self.ctx.rectangle(cx - r, cy - r, 2 * r, 2 * r)
                self.ctx.fill() if pch == 15 else self.ctx.stroke()
            elif pch == 3:  # plus
                self.ctx.set_line_width(1.0)
                self.ctx.move_to(cx - r, cy); self.ctx.line_to(cx + r, cy)
                self.ctx.move_to(cx, cy - r); self.ctx.line_to(cx, cy + r)
                self.ctx.stroke()
            elif pch == 4:  # x
                self.ctx.set_line_width(1.0)
                self.ctx.move_to(cx - r, cy - r); self.ctx.line_to(cx + r, cy + r)
                self.ctx.move_to(cx - r, cy + r); self.ctx.line_to(cx + r, cy - r)
                self.ctx.stroke()
            else:  # default solid dot
                self.ctx.arc(cx, cy, max(0.5, r * 0.5), 0, 2 * math.pi)
                self.ctx.fill()

    def text(self, x, y, label, *, color="black", size: float = 12.0,
             halign: str = "center", valign: str = "middle", angle: float = 0.0,
             family: str = "DejaVu Sans"):
        xs = np.atleast_1d(self._map_x(x))
        ys = np.atleast_1d(self._map_y(y))
        labels = label if isinstance(label, (list, np.ndarray)) else [label] * len(xs)
        self.ctx.select_font_face(family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.ctx.set_font_size(size)
        for i in range(len(xs)):
            col = color[i] if isinstance(color, (list, np.ndarray)) else color
            self._set_color(col)
            s = str(labels[i])
            xb, yb, w, h, _, _ = self.ctx.text_extents(s)
            dx = -xb if halign == "left" else (-w / 2 - xb if halign == "center" else -w - xb)
            dy = (-yb if valign == "bottom"
                  else (-yb - h / 2 if valign == "middle"
                        else (-yb - h if valign == "top" else 0)))
            cx, cy = float(xs[i]), float(ys[i])
            self.ctx.save()
            self.ctx.translate(cx, cy)
            if angle:
                self.ctx.rotate(-math.radians(angle))
            self.ctx.move_to(dx, dy)
            self.ctx.show_text(s)
            self.ctx.restore()

    def text_extents(self, label: str, size: float = 12.0,
                     family: str = "DejaVu Sans") -> tuple[float, float]:
        """Return text (width, height) in **data y units**."""
        self.ctx.select_font_face(family, cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.ctx.set_font_size(size)
        _, _, w, h, _, _ = self.ctx.text_extents(label)
        # w is in pixels; convert to data units
        wx = w / self.width_px * (self.xlim[1] - self.xlim[0])
        hy = h / self.height_px * (self.ylim[1] - self.ylim[0])
        return wx, hy

    def polygon(self, x, y, *, fill=None, border="black", lwd: float = 1.0):
        xs = self._map_x(x)
        ys = self._map_y(y)
        if len(xs) < 2:
            return
        self.ctx.move_to(float(xs[0]), float(ys[0]))
        for i in range(1, len(xs)):
            self.ctx.line_to(float(xs[i]), float(ys[i]))
        self.ctx.close_path()
        if fill is not None:
            self._set_color(fill)
            self.ctx.fill_preserve()
        if border is not None:
            self._set_color(border)
            self.ctx.set_line_width(lwd)
            self.ctx.stroke()
        else:
            self.ctx.new_path()

    def arrows(self, x0, y0, x1, y1, *, color="black", lwd: float = 1.0,
               head_length: float = 6.0, head_angle: float = 30.0):
        x0a = np.atleast_1d(self._map_x(x0))
        x1a = np.atleast_1d(self._map_x(x1))
        y0a = np.atleast_1d(self._map_y(y0))
        y1a = np.atleast_1d(self._map_y(y1))
        n = max(len(x0a), len(x1a), len(y0a), len(y1a))
        self.ctx.set_line_width(lwd)
        for i in range(n):
            col = color[i] if isinstance(color, (list, np.ndarray)) else color
            self._set_color(col)
            ax, ay = float(x0a[i % len(x0a)]), float(y0a[i % len(y0a)])
            bx, by = float(x1a[i % len(x1a)]), float(y1a[i % len(y1a)])
            self.ctx.move_to(ax, ay); self.ctx.line_to(bx, by)
            self.ctx.stroke()
            theta = math.atan2(by - ay, bx - ax)
            ang = math.radians(head_angle)
            for sign in (1, -1):
                hx = bx - head_length * math.cos(theta + sign * ang)
                hy = by - head_length * math.sin(theta + sign * ang)
                self.ctx.move_to(bx, by); self.ctx.line_to(hx, hy)
                self.ctx.stroke()

    def _set_dash(self, lty):
        if lty in (None, "solid", 1, "1"):
            self.ctx.set_dash([])
        elif lty in ("dashed", 2):
            self.ctx.set_dash([6, 4])
        elif lty in ("dotted", 3):
            self.ctx.set_dash([1, 3])
        elif lty in ("dotdash", 4):
            self.ctx.set_dash([1, 3, 6, 3])
        else:
            self.ctx.set_dash([])

    # ------------------------------------------------------------------
    # Clipping
    # ------------------------------------------------------------------
    def push_clip(self, xleft, ybottom, xright, ytop):
        xl = float(self._map_x(xleft))
        xr = float(self._map_x(xright))
        yt = float(self._map_y(ytop))
        yb = float(self._map_y(ybottom))
        self.ctx.save()
        self.ctx.rectangle(xl, yt, xr - xl, yb - yt)
        self.ctx.clip()

    def pop_clip(self) -> None:
        self.ctx.restore()

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    def save(self, path: str) -> None:
        if self.surface_kind == "png":
            self.surface.write_to_png(path)
        else:
            self.surface.finish()
            with open(path, "wb") as f:
                f.write(self._buffer.getvalue())

    def to_bytes(self) -> bytes:
        if self.surface_kind == "png":
            buf = io.BytesIO()
            self.surface.write_to_png(buf)
            return buf.getvalue()
        self.surface.finish()
        return self._buffer.getvalue()
