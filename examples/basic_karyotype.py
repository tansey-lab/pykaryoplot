"""Reproduce a few of karyoploteR's tutorial figures."""
from pathlib import Path

import numpy as np

import pykaryoplot as pk


OUT = Path(__file__).parent / "out"
OUT.mkdir(exist_ok=True)


def example_1_simplest():
    kp = pk.plot_karyotype("hg19", main="hg19 karyotype")
    kp.save(str(OUT / "01_simple.png"))


def example_2_two_panels_with_points():
    rng = np.random.default_rng(1000)
    n = 4000
    chroms = rng.choice(["chr1", "chr2"], size=n)
    starts = rng.integers(1, 240_000_000, size=n)
    ys = rng.normal(0.5, 0.1, size=n).clip(0, 1)

    kp = pk.plot_karyotype("hg19", plot_type=2,
                           chromosomes=["chr1", "chr2"])
    pk.kp_data_background(kp, color="#E8F5E8", data_panel=1)
    pk.kp_data_background(kp, color="#E8E8F5", data_panel=2)
    pk.kp_points(kp, chrom=chroms, x=starts, y=ys, data_panel=1, color="#226622", size=1.5)
    pk.kp_points(kp, chrom=chroms, x=starts, y=ys, data_panel=2, color="#222266", size=1.5)
    pk.kp_axis(kp, data_panel=1, n_ticks=3)
    pk.kp_axis(kp, data_panel=2, n_ticks=3)
    kp.save(str(OUT / "02_two_panels.png"))


def example_3_density():
    rng = np.random.default_rng(7)
    starts = rng.integers(1, 240_000_000, size=10_000)
    regions = list(zip(rng.choice(["chr1", "chr2", "chr3"], size=10_000), starts, starts + 1_000))
    kp = pk.plot_karyotype("hg19", plot_type=4,
                           chromosomes=["chr1", "chr2", "chr3"])
    pk.kp_data_background(kp, color="#FAFAFA")
    pk.kp_plot_density(kp, regions, window_size=2_000_000)
    pk.kp_add_chromosome_separators(kp)
    kp.save(str(OUT / "03_density.png"))


def example_4_manhattan():
    rng = np.random.default_rng(123)
    n = 5000
    chroms = rng.choice(["chr1", "chr2", "chr3", "chr4"], size=n)
    pos = rng.integers(1, 200_000_000, size=n)
    pvals = rng.beta(0.5, 100, size=n)
    pvals[rng.choice(n, 8, replace=False)] = 10 ** -rng.uniform(7, 10, size=8)
    kp = pk.plot_karyotype("hg19", plot_type=4,
                           chromosomes=["chr1", "chr2", "chr3", "chr4"])
    pk.kp_data_background(kp, color="#FAFAFA")
    pk.kp_plot_manhattan(kp, chrom=chroms, x=pos, pval=pvals)
    pk.kp_axis(kp, ymin=0, ymax=10, n_ticks=5)
    kp.save(str(OUT / "04_manhattan.png"))


def example_5_markers():
    kp = pk.plot_karyotype("hg19", plot_type=1, chromosomes=["chr17"])
    markers = [("chr17", 41_276_135, 41_277_500),  # BRCA1-ish
               ("chr17",  7_661_779,  7_687_550),  # TP53
               ("chr17", 37_844_393, 37_886_679)]  # ERBB2
    names = ["BRCA1", "TP53", "ERBB2"]
    pk.kp_plot_markers(kp, chrom=[m[0] for m in markers],
                       x=[(m[1] + m[2]) / 2 for m in markers],
                       labels=names, y=0.6)
    kp.save(str(OUT / "05_markers.png"))


if __name__ == "__main__":
    example_1_simplest()
    example_2_two_panels_with_points()
    example_3_density()
    example_4_manhattan()
    example_5_markers()
    print(f"Wrote 5 PNGs to {OUT}")
