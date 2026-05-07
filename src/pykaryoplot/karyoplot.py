"""Top-level :class:`KaryoPlot` object plus :func:`plot_karyotype`."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Optional, Sequence

import numpy as np

from .backend_cairo import CairoBackend
from .coord_change import CoordChange, build_coord_change
from .cytobands import filter_canonical, get_cytobands, get_genome
from .plot_params import get_default_plot_params
from .ranges import GRanges, to_granges


@dataclass
class KaryoPlot:
    plot_params: dict
    plot_type: int
    genome: GRanges
    chromosomes: list[str]
    chromosome_lengths: dict[str, int]
    chromosome_starts: dict[str, int]
    cytobands: Optional[GRanges]
    plot_region: GRanges
    zoom: bool
    coord: CoordChange
    backend: CairoBackend
    genome_name: str = "custom"

    @property
    def coord_change(self) -> Callable:
        return self.coord.coord_change

    @property
    def ideogram_mid(self) -> Callable:
        return self.coord.ideogram_mid

    @property
    def chromosome_height(self) -> float:
        return self.coord.chr_height

    def save(self, path: str) -> None:
        self.backend.save(path)


def _resolve_genome(genome) -> tuple[GRanges, str]:
    if isinstance(genome, GRanges):
        return genome, "custom"
    if isinstance(genome, str):
        gr = get_genome(genome)
        if gr is None:
            raise ValueError(f"Unknown genome {genome!r}; bundled genomes: "
                             "hg19, hg38, mm9, mm10, mm39, dm6, sacCer3, danRer11, rn6. "
                             "Pass a GRanges directly for custom genomes.")
        return gr, genome
    return to_granges(genome), "custom"


def plot_karyotype(genome="hg19", *, plot_type: int = 1,
                   chromosomes: Sequence[str] | str | None = "auto",
                   zoom: Any = None,
                   cytobands: Any = None,
                   plot_params: dict | None = None,
                   main: str | None = None,
                   width_px: int = 1200,
                   height_px: int = 800,
                   surface: str = "png",
                   data2_invert: bool = True,
                   ideogram_plotter: Callable | None = None,
                   labels_plotter: Callable | None = None) -> KaryoPlot:
    """Create a karyotype plot. Mirrors ``plotKaryotype`` from karyoploteR.

    See ``R/plotKaryotype.R`` for full parameter semantics. ``ideogram_plotter``
    and ``labels_plotter`` default to :func:`kp_add_cytobands` and
    :func:`kp_add_chromosome_names`; pass ``None`` to skip either step.
    """
    if plot_params is None:
        plot_params = get_default_plot_params(plot_type)
    plot_params.setdefault("data2_invert", True)
    if not data2_invert:
        plot_params["data2_invert"] = False

    gr_genome, genome_name = _resolve_genome(genome)

    # Zoom takes precedence and constrains to one chromosome
    zoom_gr = None
    if zoom is not None:
        zoom_gr = to_granges(zoom)
        if len(zoom_gr) > 1:
            zoom_gr = zoom_gr[[0]]
        chromosomes = [str(zoom_gr.chrom[0])]

    # Chromosome filtering
    if chromosomes is not None and chromosomes != "all":
        if isinstance(chromosomes, str) and chromosomes in ("auto", "canonical"):
            gr_genome = filter_canonical(gr_genome)
        elif isinstance(chromosomes, str) and chromosomes == "autosomal":
            keep = [c for c in gr_genome.chrom
                    if str(c).lower() not in ("chrx", "chry", "x", "y", "chrm", "mt")]
            gr_genome = filter_canonical(gr_genome).filter_chromosomes(keep)
        else:
            keep = list(chromosomes)
            gr_genome = gr_genome.filter_chromosomes(keep)
            # preserve user-specified order
            order = {c: i for i, c in enumerate(keep)}
            idx = np.argsort([order.get(c, 1e9) for c in gr_genome.chrom])
            gr_genome = gr_genome[idx]

    if len(gr_genome) == 0:
        raise ValueError("No chromosomes left after filtering.")

    # Cytobands: explicit > genome-name lookup > None
    cyto = None
    if cytobands is not None:
        cyto = to_granges(cytobands)
    elif isinstance(genome, str):
        cyto = get_cytobands(genome)
    if cyto is not None:
        cyto = cyto.filter_chromosomes(list(gr_genome.chrom))

    chrom_list = [str(c) for c in gr_genome.chrom]
    chr_starts = {c: int(s) for c, s in zip(chrom_list, gr_genome.start)}
    chr_ends = {c: int(e) for c, e in zip(chrom_list, gr_genome.end)}

    plot_region = zoom_gr if zoom_gr is not None else gr_genome

    coord = build_coord_change(plot_type, chrom_list, chr_starts, chr_ends, plot_params)

    backend = CairoBackend.create(coord.plot_xlim, coord.plot_ylim,
                                  width_px=width_px, height_px=height_px, kind=surface)

    kp = KaryoPlot(
        plot_params=plot_params,
        plot_type=plot_type,
        genome=gr_genome,
        chromosomes=chrom_list,
        chromosome_lengths={c: chr_ends[c] - chr_starts[c] for c in chrom_list},
        chromosome_starts=chr_starts,
        cytobands=cyto,
        plot_region=plot_region,
        zoom=zoom_gr is not None,
        coord=coord,
        backend=backend,
        genome_name=genome_name,
    )

    # Default plotters — imported lazily to avoid cycles
    from .ideogram import kp_add_cytobands, kp_add_chromosome_names, kp_add_main_title

    if ideogram_plotter is None:
        ideogram_plotter = kp_add_cytobands
    if labels_plotter is None:
        labels_plotter = kp_add_chromosome_names

    if ideogram_plotter is not None:
        ideogram_plotter(kp)
    if labels_plotter is not None:
        labels_plotter(kp)
    if main is not None:
        kp_add_main_title(kp, main)

    return kp
