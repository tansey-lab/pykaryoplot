"""Lightweight GRanges-like container plus coercion helpers."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping, Sequence

import numpy as np


@dataclass
class GRanges:
    """A minimal genomic-ranges container.

    Coordinates are 1-based, inclusive on both ends, matching karyoploteR /
    Bioconductor conventions. ``mcols`` holds optional per-range metadata
    columns (e.g. ``y``, ``value``, ``name``, ``gie_stain``).
    """

    chrom: np.ndarray
    start: np.ndarray
    end: np.ndarray
    mcols: dict[str, np.ndarray] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.chrom = np.asarray(self.chrom, dtype=object)
        self.start = np.asarray(self.start, dtype=np.int64)
        self.end = np.asarray(self.end, dtype=np.int64)
        n = len(self.chrom)
        if len(self.start) != n or len(self.end) != n:
            raise ValueError("chrom/start/end length mismatch")
        for k, v in list(self.mcols.items()):
            arr = np.asarray(v)
            if len(arr) != n:
                raise ValueError(f"mcol {k!r} length mismatch")
            self.mcols[k] = arr

    def __len__(self) -> int:
        return len(self.chrom)

    @property
    def midpoints(self) -> np.ndarray:
        return (self.start + self.end) / 2.0

    def __getitem__(self, key) -> "GRanges":
        if isinstance(key, str):
            return self.mcols[key]
        idx = key
        if isinstance(idx, np.ndarray) and idx.dtype == bool:
            idx = np.where(idx)[0]
        return GRanges(
            chrom=self.chrom[idx],
            start=self.start[idx],
            end=self.end[idx],
            mcols={k: v[idx] for k, v in self.mcols.items()},
        )

    def filter_chromosomes(self, keep: Sequence[str]) -> "GRanges":
        keep_set = set(keep)
        mask = np.array([c in keep_set for c in self.chrom], dtype=bool)
        return self[mask]


def to_granges(obj: Any) -> GRanges:
    """Coerce common inputs into a :class:`GRanges`.

    Accepts: an existing GRanges, a dict-like with chrom/start/end columns,
    a pandas DataFrame, a pyranges.PyRanges, or a sequence of (chrom, start,
    end[, value]) tuples.
    """
    if isinstance(obj, GRanges):
        return obj

    # pandas / pyranges duck-typed: has columns
    cols = getattr(obj, "columns", None)
    if cols is not None:
        # pyranges has Chromosome/Start/End; pandas typically chrom/start/end
        col_map = {c.lower(): c for c in cols}
        chrom_col = col_map.get("chrom") or col_map.get("chromosome") or col_map.get("seqnames")
        start_col = col_map.get("start")
        end_col = col_map.get("end")
        if chrom_col is None or start_col is None or end_col is None:
            raise ValueError(f"DataFrame must have chrom/start/end columns; got {list(cols)}")
        df = obj.df if hasattr(obj, "df") else obj  # pyranges → DataFrame
        mcols = {}
        for c in cols:
            if c in (chrom_col, start_col, end_col):
                continue
            mcols[c.lower() if c.lower() in ("y", "value", "name", "gie_stain", "giestain") else c] = np.asarray(df[c])
        # rename giestain -> gie_stain
        if "giestain" in mcols and "gie_stain" not in mcols:
            mcols["gie_stain"] = mcols.pop("giestain")
        return GRanges(np.asarray(df[chrom_col]),
                       np.asarray(df[start_col]),
                       np.asarray(df[end_col]),
                       mcols)

    if isinstance(obj, Mapping):
        return GRanges(
            obj["chrom"],
            obj["start"],
            obj["end"],
            {k: v for k, v in obj.items() if k not in ("chrom", "start", "end")},
        )

    if isinstance(obj, Iterable):
        rows = list(obj)
        if not rows:
            return GRanges(np.array([], dtype=object), np.array([], dtype=np.int64), np.array([], dtype=np.int64))
        widths = {len(r) for r in rows}
        if widths == {3}:
            chrom, start, end = zip(*rows)
            return GRanges(np.array(chrom, dtype=object), np.array(start), np.array(end))
        if widths == {4}:
            chrom, start, end, val = zip(*rows)
            return GRanges(np.array(chrom, dtype=object), np.array(start), np.array(end),
                           {"value": np.array(val)})
        raise ValueError("Tuple input must be (chrom,start,end) or (chrom,start,end,value)")

    raise TypeError(f"Cannot coerce {type(obj).__name__} to GRanges")
