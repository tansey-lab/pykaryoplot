import os
from pathlib import Path

import numpy as np

import pykaryoplot as pk


def test_render_basic_hg19(tmp_path):
    kp = pk.plot_karyotype("hg19", chromosomes=["chr1", "chr2", "chrX"])
    out = tmp_path / "basic.png"
    kp.save(str(out))
    assert out.exists() and out.stat().st_size > 1000


def test_render_with_data(tmp_path):
    rng = np.random.default_rng(0)
    n = 500
    chroms = rng.choice(["chr1", "chr2"], size=n)
    starts = rng.integers(1, 100_000_000, size=n)
    ys = rng.uniform(0, 1, size=n)
    kp = pk.plot_karyotype("hg19", plot_type=2,
                           chromosomes=["chr1", "chr2"])
    pk.kp_data_background(kp, color="#F0F0F0", data_panel=1)
    pk.kp_data_background(kp, color="#F0F0F8", data_panel=2)
    pk.kp_points(kp, chrom=chroms, x=starts, y=ys, data_panel=1, size=2)
    pk.kp_lines(kp, chrom=chroms, x=starts, y=ys, data_panel=2, color="#3333AA")
    out = tmp_path / "data.png"
    kp.save(str(out))
    assert out.exists() and out.stat().st_size > 1000


def test_single_line_plot_type_4(tmp_path):
    kp = pk.plot_karyotype("hg19", plot_type=4,
                           chromosomes=["chr1", "chr2", "chr3"])
    pk.kp_data_background(kp, color="#EEEEEE")
    out = tmp_path / "type4.png"
    kp.save(str(out))
    assert out.exists()


def test_density_plot(tmp_path):
    rng = np.random.default_rng(42)
    starts = rng.integers(1, 240_000_000, size=2000)
    regions = list(zip(["chr1"] * 2000, starts, starts + 10_000))
    kp = pk.plot_karyotype("hg19", chromosomes=["chr1"])
    pk.kp_data_background(kp, color="#FAFAFA")
    pk.kp_plot_density(kp, regions, window_size=2_000_000)
    out = tmp_path / "density.png"
    kp.save(str(out))
    assert out.exists() and out.stat().st_size > 1000
