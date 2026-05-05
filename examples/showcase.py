"""Showcase plots demonstrating more pykaryoplot features."""
from pathlib import Path

import numpy as np

import pykaryoplot as pk

OUT = Path(__file__).parent / "out"
OUT.mkdir(exist_ok=True)


def panel_with_tracks():
    """Stacked autotracks: 4 separate signal tracks in a single data panel."""
    rng = np.random.default_rng(2026)
    kp = pk.plot_karyotype("hg38", plot_type=1,
                           chromosomes=["chr1", "chr2", "chr3"],
                           main="Stacked autotracks (hg38)",
                           width_px=1400, height_px=600)
    pk.kp_data_background(kp, color="#FAFAFA")
    palette = ["#1F77B4", "#FF7F0E", "#2CA02C", "#D62728"]
    for i in range(1, 5):
        r0, r1 = pk.prepare_params.autotrack(i, 4)  # type: ignore[attr-defined]
        n = 300
        chroms = rng.choice(["chr1", "chr2", "chr3"], size=n)
        starts = rng.integers(1, 240_000_000, size=n)
        ys = rng.uniform(0, 1, size=n)
        pk.kp_points(kp, chrom=chroms, x=starts, y=ys, r0=r0, r1=r1,
                     color=palette[i - 1], size=1.5)
        pk.kp_add_color_rect(kp, color=palette[i - 1], r0=r0, r1=r1)
        pk.kp_add_labels(kp, f"track {i}", r0=r0, r1=r1, cex=0.7)
    kp.save(str(OUT / "10_tracks.png"))


def heatmap_demo():
    """Sliding-window heatmap with col_by_value across hg38."""
    rng = np.random.default_rng(0xC0FFEE)
    kp = pk.plot_karyotype("hg38", plot_type=4,
                           chromosomes=[f"chr{i}" for i in range(1, 6)],
                           main="GC-content-like heatmap",
                           width_px=1400, height_px=400)
    chr_arr, x0_arr, x1_arr, y_arr = [], [], [], []
    for c in kp.chromosomes:
        s0 = kp.chromosome_starts[c]
        e0 = s0 + kp.chromosome_lengths[c]
        win = 500_000
        starts = np.arange(s0, e0, win)
        ends = np.minimum(starts + win, e0)
        # smoothed random signal
        sig = np.cumsum(rng.normal(0, 1, size=len(starts)))
        sig = (sig - sig.min()) / (sig.max() - sig.min() + 1e-9)
        chr_arr.extend([c] * len(starts))
        x0_arr.extend(starts.tolist())
        x1_arr.extend(ends.tolist())
        y_arr.extend(sig.tolist())
    pk.kp_heatmap(kp, chrom=np.array(chr_arr, dtype=object),
                  x0=np.array(x0_arr), x1=np.array(x1_arr),
                  y=np.array(y_arr),
                  colors=("#08306B", "#F7FBFF", "#A50F15"))
    pk.kp_add_chromosome_separators(kp, color="#444444")
    kp.save(str(OUT / "11_heatmap.png"))


def coverage_and_regions():
    """kp_plot_coverage on top, kp_plot_regions stacked underneath."""
    rng = np.random.default_rng(99)
    n = 600
    chroms = rng.choice(["chr19", "chr20", "chr21", "chr22"], size=n)
    starts = rng.integers(1, 60_000_000, size=n)
    widths = rng.integers(100_000, 5_000_000, size=n)
    regions = list(zip(chroms, starts, starts + widths))

    kp = pk.plot_karyotype("hg19", plot_type=2,
                           chromosomes=["chr19", "chr20", "chr21", "chr22"],
                           main="Coverage vs region tracks",
                           width_px=1200, height_px=500)
    pk.kp_data_background(kp, color="#F5FAF5", data_panel=1)
    pk.kp_data_background(kp, color="#F5F5FA", data_panel=2)
    pk.kp_plot_coverage(kp, regions, data_panel=1, color="#226622")
    pk.kp_plot_regions(kp, regions, data_panel=2, color="#5577CC", border="#223388")
    pk.kp_add_labels(kp, "coverage", data_panel=1, cex=0.8)
    pk.kp_add_labels(kp, "regions",  data_panel=2, cex=0.8)
    kp.save(str(OUT / "12_coverage_regions.png"))


def ribbon_and_lines():
    """Confidence-band ribbon + center line on a zoomed plot."""
    rng = np.random.default_rng(3)
    kp = pk.plot_karyotype("hg38", plot_type=1, chromosomes=["chr7"],
                           main="Confidence ribbon (chr7)",
                           width_px=1400, height_px=350)
    pk.kp_data_background(kp, color="#FAFAFA")
    n = 400
    x = np.linspace(1, kp.chromosome_lengths["chr7"], n).astype(np.int64)
    base = 0.5 + 0.3 * np.sin(x / 1e7)
    noise = rng.normal(0, 0.05, n)
    center = base + noise
    lo = center - 0.1
    hi = center + 0.1
    chrom = np.array(["chr7"] * n, dtype=object)
    pk.kp_plot_ribbon(kp, chrom=chrom, x=x, y_lo=lo, y_hi=hi,
                      ymin=0, ymax=1, color="#88AAEE66", border="#3355BB")
    pk.kp_lines(kp, chrom=chrom, x=x, y=center, ymin=0, ymax=1,
                color="#1A2B66", lwd=1.5)
    pk.kp_axis(kp, ymin=0, ymax=1, n_ticks=5)
    kp.save(str(OUT / "13_ribbon.png"))


def yeast_full_genome():
    kp = pk.plot_karyotype("sacCer3", plot_type=4,
                           main="S. cerevisiae genome",
                           width_px=1400, height_px=300)
    pk.kp_add_chromosome_separators(kp, color="#888888")
    kp.save(str(OUT / "14_yeast.png"))


if __name__ == "__main__":
    panel_with_tracks()
    heatmap_demo()
    coverage_and_regions()
    ribbon_and_lines()
    yeast_full_genome()
    print(f"Wrote 5 showcase PNGs to {OUT}")
