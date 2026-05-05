"""Bundled chromosome length and cytoband tables (sourced from UCSC)."""
from __future__ import annotations

import csv
from importlib import resources
from typing import Optional

import numpy as np

from .ranges import GRanges


def _data_files(subdir: str) -> dict[str, "resources.abc.Traversable"]:
    root = resources.files("pykaryoplot").joinpath("data", subdir)
    return {p.stem: p for p in root.iterdir() if p.suffix == ".tsv"}


def _read_tsv(path) -> list[dict]:
    with path.open("r") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        return list(rdr)


def get_genome(name: str) -> Optional[GRanges]:
    """Return a GRanges of chromosome ranges for a UCSC-style genome name."""
    files = _data_files("genomes")
    if name not in files:
        return None
    rows = _read_tsv(files[name])
    chrom = np.array([r["chrom"] for r in rows], dtype=object)
    start = np.array([int(r["start"]) for r in rows], dtype=np.int64)
    end = np.array([int(r["end"]) for r in rows], dtype=np.int64)
    return GRanges(chrom, start, end)


def get_cytobands(name: str) -> Optional[GRanges]:
    files = _data_files("cytobands")
    if name not in files:
        return None
    rows = _read_tsv(files[name])
    if not rows:
        return None
    chrom = np.array([r["chrom"] for r in rows], dtype=object)
    start = np.array([int(r["start"]) + 1 for r in rows], dtype=np.int64)  # UCSC is 0-based
    end = np.array([int(r["end"]) for r in rows], dtype=np.int64)
    name_arr = np.array([r["name"] for r in rows], dtype=object)
    stain = np.array([r["gie_stain"] for r in rows], dtype=object)
    return GRanges(chrom, start, end, {"name": name_arr, "gie_stain": stain})


def filter_canonical(genome: GRanges) -> GRanges:
    """Heuristic: drop chrUn*, *_random, *_alt, *_fix, *_hap*, chrM."""
    keep = []
    for c in genome.chrom:
        cl = str(c).lower()
        if cl.startswith("chrun") or cl.startswith("un_"):
            continue
        if any(s in cl for s in ("_random", "_alt", "_fix", "_hap", "_ctg")):
            continue
        if cl in ("chrm", "mt", "chrmt"):
            continue
        keep.append(c)
    out = genome.filter_chromosomes(keep)
    # natural sort: numeric chromosomes first in numeric order, then X/Y/M
    def sort_key(c):
        s = str(c).lower().removeprefix("chr")
        try:
            return (0, int(s))
        except ValueError:
            return (1, {"x": 23, "y": 24, "m": 25, "mt": 25}.get(s, 99))
    keys = [sort_key(c) for c in out.chrom]
    idx = sorted(range(len(keys)), key=lambda i: keys[i])
    return out[np.array(idx)]
