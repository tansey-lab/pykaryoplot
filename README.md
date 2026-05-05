# pykaryoplot

A Python port of the R/Bioconductor package
[**karyoploteR**](https://github.com/bernatgel/karyoploteR), rendered with
[pycairo](https://pycairo.readthedocs.io/). Build customizable linear-genome
karyotype plots and overlay arbitrary genomic data — points, lines, bars,
heatmaps, ribbons, density tracks, Manhattan plots, and more — without
leaving the Python ecosystem.

Pure Python at runtime: depends only on `pycairo` and `numpy`. No R, no
Bioconductor, no system bioinformatics toolchain required for the core
plotting surface.

## Install

```bash
pip install pykaryoplot
# or, for the optional input-coercion / IO extras:
pip install "pykaryoplot[pyranges,pandas]"
```

Bundled UCSC chromosome and cytoband tables for `hg19`, `hg38`, `mm9`,
`mm10`, `mm39`, `dm6`, `sacCer3`, `danRer11`, and `rn6`. Custom genomes
are supported by passing your own chromosome lengths and cytoband table.

## Quickstart

```python
import numpy as np
import pykaryoplot as pk

rng = np.random.default_rng(0)
n = 4000
chroms = rng.choice(["chr1", "chr2"], size=n)
pos    = rng.integers(1, 240_000_000, size=n)
y      = rng.normal(0.5, 0.1, size=n).clip(0, 1)

kp = pk.plot_karyotype("hg19", plot_type=2, chromosomes=["chr1", "chr2"])
pk.kp_data_background(kp, color="#E8F5E8", data_panel=1)
pk.kp_data_background(kp, color="#E8E8F5", data_panel=2)
pk.kp_points(kp, chrom=chroms, x=pos, y=y, data_panel=1, color="#226622")
pk.kp_lines (kp, chrom=chroms, x=pos, y=y, data_panel=2, color="#222266")
pk.kp_axis  (kp, data_panel=1, n_ticks=3)
pk.kp_axis  (kp, data_panel=2, n_ticks=3)
kp.save("two_panels.png")
```

See `examples/basic_karyotype.py` and `examples/showcase.py` for richer
demos: stacked autotracks, GC-style heatmaps, coverage + region tracks,
confidence ribbons, density / rainfall / Manhattan plots, marker labels.

## API

`plot_karyotype(genome, plot_type, chromosomes, zoom, cytobands, plot_params, main, ...)`
returns a `KaryoPlot` you pass to every `kp_*` call. The same seven
`plot.type` layouts as karyoploteR are supported (1, 2, 6 stack each
chromosome on its own row; 3, 4, 5, 7 lay all chromosomes in a single row).

### Implemented

- **Ideogram / annotation** — `kp_add_cytobands`,
  `kp_add_cytobands_as_line`, `kp_add_cytoband_labels`,
  `kp_add_chromosome_names`, `kp_add_chromosome_separators`,
  `kp_add_base_numbers`, `kp_add_main_title`, `kp_add_labels`,
  `kp_add_color_rect`
- **Data primitives** — `kp_points`, `kp_lines`, `kp_segments`,
  `kp_arrows`, `kp_text`, `kp_abline`, `kp_rect`, `kp_polygon`,
  `kp_bars`, `kp_area`, `kp_heatmap`, `kp_data_background`, `kp_axis`
- **High-level plotters** — `kp_plot_regions`, `kp_plot_density`,
  `kp_plot_coverage`, `kp_plot_markers`, `kp_plot_2lines`,
  `kp_plot_rainfall`, `kp_plot_manhattan`, `kp_plot_ribbon`
- **Color helpers** — `darker`, `lighter`, `transparent`,
  `col_by_chr`, `col_by_value`, plus the three karyoploteR cytoband
  schemas (`circos`, `biovizbase`, `only.centromeres`)

### Deferred (planned for v2)

- BAM-backed plotters (`kp_plot_bam_coverage`, `kp_plot_bam_density`) —
  needs `pysam`
- BigWig-backed plotters (`kp_plot_bigwig`) — needs `pyBigWig`
- Gene / transcript plotters (`kp_plot_genes`, `kp_plot_transcripts`)
- LOESS smoothing (`kp_plot_loess`), horizon plots (`kp_plot_horizon`),
  full Bezier `kp_plot_links`

## Inputs

`kp_*` functions accept a few input shapes for genomic data:

- the bundled `pykaryoplot.GRanges` (chrom, start, end, mcols dict),
- a `pyranges.PyRanges` (column names auto-detected),
- a `pandas.DataFrame` with `chrom`, `start`, `end` columns,
- a list of `(chrom, start, end[, value])` tuples,
- explicit `chrom=`, `x=`, `y=` arrays.

## Output

PNG, SVG, or PDF (`plot_karyotype(..., surface="svg")`); `kp.save(path)`.

## Acknowledgments

This package is a port of [karyoploteR](https://github.com/bernatgel/karyoploteR)
by Bernat Gel. Bundled cytoband data is from the
[UCSC Genome Browser](https://genome.ucsc.edu/).

## License

[Artistic-2.0](LICENSE), matching upstream `karyoploteR`.
