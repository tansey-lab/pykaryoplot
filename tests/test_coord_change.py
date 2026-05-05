import numpy as np

from pykaryoplot.coord_change import build_coord_change
from pykaryoplot.plot_params import get_default_plot_params


def _two_chr_coord(plot_type=1):
    pp = get_default_plot_params(plot_type)
    chr_names = ["chr1", "chr2"]
    chr_starts = {"chr1": 1, "chr2": 1}
    chr_ends = {"chr1": 100, "chr2": 200}
    return build_coord_change(plot_type, chr_names, chr_starts, chr_ends, pp), pp


def test_x_at_chromosome_start_is_left_margin():
    cc, pp = _two_chr_coord(1)
    xp, _ = cc.coord_change(chrom=["chr1"], x=[1])
    assert abs(float(xp[0]) - pp["leftmargin"]) < 1e-12


def test_x_proportional_to_chromosome_length():
    cc, pp = _two_chr_coord(1)
    # chr2 is 200bp; midpoint x=100 should be at half its visible width (chr2 is the longest, sets max_chr_len)
    xp, _ = cc.coord_change(chrom=["chr2"], x=[100])
    expected = pp["leftmargin"] + (99 / 199) * (1 - pp["leftmargin"] - pp["rightmargin"])
    assert abs(float(xp[0]) - expected) < 1e-9


def test_single_line_layout_lays_chromosomes_in_a_row():
    cc, pp = _two_chr_coord(3)
    xp_a, _ = cc.coord_change(chrom=["chr1"], x=[100])
    xp_b, _ = cc.coord_change(chrom=["chr2"], x=[1])
    assert float(xp_a[0]) < float(xp_b[0])


def test_y_in_data_panel_1_is_above_ideogram():
    cc, pp = _two_chr_coord(2)
    _, y_ideo = cc.coord_change(chrom=["chr1"], y=[0.5], data_panel="ideogram")
    _, y_d1 = cc.coord_change(chrom=["chr1"], y=[0.5], data_panel=1)
    _, y_d2 = cc.coord_change(chrom=["chr1"], y=[0.5], data_panel=2)
    assert float(y_d2[0]) < float(y_ideo[0]) < float(y_d1[0])
