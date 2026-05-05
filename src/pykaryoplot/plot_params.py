"""Default plot.params for the seven karyoploteR plot types.

Verbatim port of ``R/getDefaultPlotParams.R``.
"""
from __future__ import annotations

VALID_PLOT_TYPES = (1, 2, 3, 4, 5, 6, 7)


def get_default_plot_params(plot_type: int) -> dict:
    if plot_type not in VALID_PLOT_TYPES:
        raise ValueError(f"plot_type must be one of {VALID_PLOT_TYPES}, got {plot_type}")

    if plot_type == 1:
        pp = dict(leftmargin=0.1, rightmargin=0.05, topmargin=120, bottommargin=100,
                  ideogramheight=50, ideogramlateralmargin=0,
                  data1height=200, data1inmargin=20, data1outmargin=20,
                  data1min=0, data1max=1,
                  data2height=0, data2inmargin=0, data2outmargin=0,
                  data2min=0, data2max=1)
    elif plot_type == 2:
        pp = dict(leftmargin=0.1, rightmargin=0.05, topmargin=120, bottommargin=100,
                  ideogramheight=50, ideogramlateralmargin=0,
                  data1height=200, data1inmargin=20, data1outmargin=20,
                  data1min=0, data1max=1,
                  data2height=200, data2inmargin=20, data2outmargin=20,
                  data2min=0, data2max=1)
    elif plot_type == 3:
        pp = dict(leftmargin=0.05, rightmargin=0.05, topmargin=30, bottommargin=30,
                  ideogramheight=10, ideogramlateralmargin=0.003,
                  data1height=200, data1inmargin=10, data1outmargin=0,
                  data1min=0, data1max=1,
                  data2height=200, data2inmargin=10, data2outmargin=0,
                  data2min=0, data2max=1)
    elif plot_type == 4:
        pp = dict(leftmargin=0.05, rightmargin=0.05, topmargin=30, bottommargin=30,
                  ideogramheight=10, ideogramlateralmargin=0.003,
                  data1height=200, data1inmargin=10, data1outmargin=0,
                  data1min=0, data1max=1,
                  data2height=0, data2inmargin=0, data2outmargin=0,
                  data2min=0, data2max=1)
    elif plot_type == 5:
        pp = dict(leftmargin=0.05, rightmargin=0.05, topmargin=30, bottommargin=30,
                  ideogramheight=10, ideogramlateralmargin=0.003,
                  data1height=0, data1inmargin=0, data1outmargin=0,
                  data1min=0, data1max=1,
                  data2height=200, data2inmargin=10, data2outmargin=0,
                  data2min=0, data2max=1)
    elif plot_type == 6:
        pp = dict(leftmargin=0.1, rightmargin=0.05, topmargin=120, bottommargin=100,
                  ideogramheight=50, ideogramlateralmargin=0,
                  data1height=0, data1inmargin=0, data1outmargin=5,
                  data1min=0, data1max=1,
                  data2height=0, data2inmargin=0, data2outmargin=5,
                  data2min=0, data2max=1)
    else:  # 7
        pp = dict(leftmargin=0.1, rightmargin=0.05, topmargin=30, bottommargin=30,
                  ideogramheight=200, ideogramlateralmargin=0.003,
                  data1height=0, data1inmargin=0, data1outmargin=5,
                  data1min=0, data1max=1,
                  data2height=0, data2inmargin=0, data2outmargin=5,
                  data2min=0, data2max=1)

    pp.update(dataideogrammin=0, dataideogrammax=1, dataallmin=0, dataallmax=1)
    return pp


SINGLE_LINE_TYPES = frozenset({3, 4, 5, 7})
STACKED_TYPES = frozenset({1, 2, 6})
