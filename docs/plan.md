# pykaryoplot — Python Port of karyoploteR

## Context

`karyoploteR` (Bernat Gel, Bioconductor) is a popular R package that draws customizable linear-genome karyotype plots and overlays arbitrary genomic data on them. We want a faithful Python port (`pykaryoplot`) so that users in the Python bioinformatics ecosystem (pysam / pyranges / pyBigWig users) can produce the same plots without leaving Python.

The R package is built on R base graphics. For the Python port we will use **pycairo** as the rendering backend, exposing PNG / SVG / PDF output via Cairo surfaces, and we will mirror karyoploteR's public API as closely as possible (snake_case'd names, e.g. `plot_karyotype`, `kp_points`).

Source package: `/Users/quinnj2/Code/karyoploteR` (50 R files, ~6k LOC).
Target package: `/Users/quinnj2/Code/pykaryoplot` (currently has only stub `pyproject.toml`, `main.py`, `README.md`).

## Architectural Mirror

The R package's public surface is a `KaryoPlot` object built by `plotKaryotype()` and a family of `kp*()` functions that take `kp` as their first argument. Internally everything routes through:

1. **plot params** — a dict of margins / heights / data-min-max for the chosen `plot.type` (1–7). See `R/getDefaultPlotParams.R`.
2. **coord-change function** — `genomic2plot(chr, x, y, data.panel)` mapping (chromosome, base-pair, value) → (plot-x, plot-y) in plot-coordinates. See `R/getCoordChangeFunctions.R`. Two underlying mappings: stacked-chromosomes (types 1, 2, 6) and single-line (types 3, 4, 5, 7).
3. **prepareParameters2 / prepareParameters4** — extract chr/x/y from GRanges *or* explicit args, recycle, scale y by `ymin/ymax/r0/r1`, filter invisible chromosomes, return filter mask. See `R/prepareParameters2.R`, `R/prepareParameters4.R`.
4. **renderer** — base-R `graphics::points/lines/rect/text/...`. We will replace this with a thin **`CairoBackend`** that exposes the same operations on a `cairo.Context`.

We will keep the same conceptual layering in Python so that future `kp*` functions can be added by writing only the genome-aware logic and delegating to `CairoBackend.points/lines/rect/...`.

## Package Layout

```
pykaryoplot/
├── pyproject.toml                       # add deps: pycairo, numpy, pyranges
├── docs/plan.md                         # this file
├── examples/
│   └── basic_karyotype.py               # smoke-test driver
├── tests/
│   ├── test_coord_change.py             # numerical equivalence to R
│   ├── test_plot_params.py
│   └── test_smoke_render.py             # render PNG, check non-empty
└── src/pykaryoplot/
    ├── __init__.py                      # re-export public API
    ├── ranges.py                        # tiny GRanges-like dataclass (chr,start,end,values)
    ├── plot_params.py                   # get_default_plot_params(plot_type)
    ├── coord_change.py                  # build_coord_change(kp) → callable
    ├── karyoplot.py                     # KaryoPlot dataclass + plot_karyotype()
    ├── backend_cairo.py                 # CairoBackend: points/lines/rect/text/polygon/clip
    ├── colors.py                        # cytoband schemas (circos/biovizbase/only.centromeres),
    │                                    # darker/lighter/transparent, col_by_chr/value/category
    ├── cytobands.py                     # bundled cytoband + chromosome-length tables
    ├── data/cytobands/                  # tsv files: hg19, hg38, mm10, mm39, dm6, ... (extracted from sysdata.rda)
    ├── data/genomes/                    # tsv files: chrom name -> length, per genome
    ├── prepare_params.py                # prepare_parameters_2 / _4, filter_params, recycle
    ├── utils.py                         # autotrack, preprocess_r0_r1, process_clipping, find_intersections
    ├── primitives.py                    # kp_points/lines/rect/segments/arrows/text/abline/bars/
    │                                    # polygon/heatmap/area/data_background/axis/add_labels/
    │                                    # add_main_title/add_color_rect
    ├── ideogram.py                      # kp_add_cytobands, kp_add_cytobands_as_line,
    │                                    # kp_add_cytoband_labels, kp_add_chromosome_names,
    │                                    # kp_add_chromosome_separators, kp_add_base_numbers
    ├── plotters.py                      # kp_plot_regions, kp_plot_density, kp_plot_coverage,
    │                                    # kp_plot_markers, kp_plot_2lines, kp_plot_rainfall,
    │                                    # kp_plot_manhattan, kp_plot_ribbon, kp_plot_loess,
    │                                    # kp_plot_horizon
    └── plotters_bio.py                  # OPTIONAL extras module (importing pysam / pyBigWig
                                         # only on demand): kp_plot_bam_coverage,
                                         # kp_plot_bam_density, kp_plot_bigwig
```

## V1 Scope (this PR)

Goal: a **working, demo-able** package that can reproduce karyoploteR's first 5 tutorial figures end-to-end on `hg19` and `hg38`.

### Must-have
- `KaryoPlot` object + `plot_karyotype(genome, plot_type, chromosomes, zoom, cytobands, plot_params, main)` for `plot_type` ∈ {1, 2, 3, 4, 5, 6, 7}.
- `get_default_plot_params(plot_type)` matching `R/getDefaultPlotParams.R` exactly.
- Coordinate-change function for both layouts (stacked & single-line), matching `R/getCoordChangeFunctions.R`.
- `CairoBackend` exposing: `points`, `lines`, `segments`, `rect`, `polygon`, `text`, `arrows`, `arc`, `clip`, `text_extents`, `save_png/svg/pdf`.
- `prepare_parameters_2` and `prepare_parameters_4` with `r0 / r1 / ymin / ymax` scaling and chromosome filtering.
- Color machinery: cytoband schemas (circos / biovizbase / only.centromeres) verbatim from `R/color.R`, plus `darker`, `lighter`, `transparent`, `col_by_chr`, `col_by_value`, `col_by_region`.
- Bundled chromosome-length and cytoband tables for `hg19`, `hg38`, `mm9`, `mm10`, `mm39`, `dm6`, `sacCer3` shipped as TSV under `src/pykaryoplot/data/`. Will be one-time-extracted from `R/sysdata.rda` (using R `saveRDS` / `Rscript -e 'd<-load("...");...'` outside the package — produces plain TSVs we ship). Custom genomes: user can pass a length dict and a cytoband DataFrame.
- Ideogram primitives: `kp_add_cytobands`, `kp_add_cytobands_as_line`, `kp_add_cytoband_labels`, `kp_add_chromosome_names`, `kp_add_chromosome_separators`, `kp_add_base_numbers`, `kp_add_main_title`, `kp_add_labels`.
- Data primitives: `kp_points`, `kp_lines`, `kp_segments`, `kp_arrows`, `kp_text`, `kp_abline`, `kp_rect`, `kp_polygon`, `kp_bars`, `kp_area`, `kp_heatmap`, `kp_data_background`, `kp_axis`, `kp_add_color_rect`.
- High-level (no heavy deps): `kp_plot_regions`, `kp_plot_density`, `kp_plot_coverage`, `kp_plot_markers`, `kp_plot_2lines`, `kp_plot_rainfall`, `kp_plot_manhattan`, `kp_plot_ribbon`.
- Public input type: any object with `.chrom / .start / .end / .name / .gie_stain / .value / .y` columns — we accept a small `GRanges`-like dataclass plus auto-coercion from `pyranges.PyRanges`, pandas DataFrame (with the standard column names), and tuples-of-(chr, start, end[, value]). Input coercion lives in `ranges.py::to_granges()`.

### Out-of-scope for v1 (explicit deferral)
- BAM-backed plotters (`kp_plot_bam_coverage`, `kp_plot_bam_density`) — needs `pysam`. Stub module `plotters_bio.py` with `NotImplementedError` and a comment pointing to v2.
- BigWig-backed plotters (`kp_plot_bigwig`) — needs `pyBigWig`. Same treatment.
- Gene/transcript plotting (`kp_plot_genes`, `kp_plot_transcripts`, `make_genes_data_from_*`) — TxDb conversion is its own project; defer to v2.
- LOESS smoothing / horizon plots (`kp_plot_loess`, `kp_plot_horizon`) — port in v2 once primitives are stable.
- Inter-chromosomal links beyond simple ribbons (`kp_plot_links` complex Bezier mode) — basic ribbon yes, full bezier defer.
- UCSC live download fallback for cytobands of unknown genomes — v1 will raise with a clear message; v2 can add a `fetch_cytobands_from_ucsc()` helper.

## Critical Files To Reference While Implementing

| Concept | R source |
|---|---|
| Public entry point & filtering | `R/plotKaryotype.R` |
| Plot params per plot.type | `R/getDefaultPlotParams.R` |
| Coordinate transforms | `R/getCoordChangeFunctions.R` (lines 14–144 dispatcher; 227–394 implementations) |
| Param normalization | `R/prepareParameters2.R`, `R/prepareParameters4.R`, `R/utils.R` |
| Cytoband rendering & colors | `R/kpAddCytobands.R`, `R/color.R` (lines 9–60), `R/getCytobandColors.R` |
| Cytoband / genome data | `R/sysdata.rda` (extract once, ship TSVs), `R/getCytobands.R` |
| Annotation primitives | `R/kpAddChromosomeNames.R`, `R/kpAddMainTitle.R`, `R/kpAddBaseNumbers.R`, `R/kpAddLabels.R`, `R/kpAxis.R`, `R/kpDataBackground.R` |
| Data primitives | `R/kpPoints.R`, `R/kpRect.R`, `R/kpLines.R`, `R/kpSegments.R`, `R/kpText.R`, `R/kpBars.R`, `R/kpArea.R`, `R/kpHeatmap.R`, `R/kpAbline.R` |
| Higher-level plotters | `R/kpPlotRegions.R`, `R/kpPlotDensity.R`, `R/kpPlotCoverage.R`, `R/kpPlotMarkers.R`, `R/kpPlotRainfall.R`, `R/kpPlotManhattan.R`, `R/kpPlot2Lines.R`, `R/kpPlotRibbon.R` |

## Implementation Order

1. **Skeleton** — `pyproject.toml` (deps: `pycairo`, `numpy`, `pyranges`, optional extras `pysam`, `pyBigWig`); empty package; `tests/` with pytest.
2. **`ranges.py`** — `GRanges` dataclass + `to_granges()` coercion.
3. **`plot_params.py`** — verbatim port of `getDefaultPlotParams`. Unit-test field-by-field against R output (one fixture per plot.type).
4. **`coord_change.py`** — port both `genomic2plot_*` paths. Unit-test on synthetic 2-chromosome genome at known input/output points (the R values are deterministic; capture them once with R, store as JSON fixtures).
5. **`backend_cairo.py`** — implement the small graphics surface. Plot-coordinate origin must match base-R: `(0,0)` lower-left, y grows upward. Cairo y-axis grows downward, so the backend flips internally.
6. **`cytobands.py`** + bundled TSVs — extract from `R/sysdata.rda` with a one-shot `scripts/extract_sysdata.R` that writes TSVs into `src/pykaryoplot/data/`.
7. **`karyoplot.py`** + `plot_karyotype()` — wire 3+4+5+6 together; render an empty plot with chromosome names. End-to-end smoke test producing a PNG.
8. **`ideogram.py`** — `kp_add_cytobands` first; visual diff against an R-produced reference PNG of `plot_karyotype("hg19")`.
9. **`prepare_params.py`** — port `prepareParameters2/4`, then port `kp_points` and `kp_rect` first as the canonical examples.
10. **`primitives.py`** — fill out the rest of the data primitives.
11. **`plotters.py`** — high-level no-dep plotters.
12. **Examples** — `examples/basic_karyotype.py` reproducing tutorial figures 1–5.

## Verification

- **Unit tests**: `pytest` covering plot params, coordinate transforms (numerical equivalence to R fixtures), color schemas, parameter normalization (recycling, r0/r1, ymin/ymax, chromosome filtering).
- **Smoke render tests**: render `plot_karyotype("hg19")` to PNG; assert file is non-empty and a known-pixel-region (a chromosome name area) is non-white.
- **Visual regression**: `examples/basic_karyotype.py` produces a PNG; we will hand-compare to a reference produced by R running the corresponding R snippet from `R/plotKaryotype.R`'s `@examples`. Acceptable tolerance: same layout & colors; sub-pixel anti-aliasing differences are fine.
- **Manual run**:
  ```bash
  uv pip install -e .
  python examples/basic_karyotype.py            # writes example1.png … example5.png
  pytest -q
  ```

## Risks / Open Questions

1. **Text metrics**: pycairo's `text_extents` measures in device units; karyoploteR uses R's `strwidth/strheight` which is in user-coordinates. We'll wrap `cairo.text_extents` with a unit-conversion helper. Worst case: small label-position discrepancies vs R; not a correctness issue.
2. **Font availability**: Cairo will fall back to a default font if not specified. We will default to `"DejaVu Sans"` to be reproducible across Linux/macOS.
3. **`sysdata.rda` extraction**: requires R installed once at packaging time, not at runtime. The script `scripts/extract_sysdata.R` (a few lines) runs only on the maintainer machine; we'll commit the produced TSVs.
4. **GRanges-style API**: Python users may pass `pyranges.PyRanges`, `pandas.DataFrame`, or our own `GRanges`. Coercion is centralized in `ranges.to_granges()`.
