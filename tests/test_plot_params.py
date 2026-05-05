import pytest

from pykaryoplot.plot_params import get_default_plot_params, VALID_PLOT_TYPES


@pytest.mark.parametrize("pt", VALID_PLOT_TYPES)
def test_each_plot_type_has_required_keys(pt):
    pp = get_default_plot_params(pt)
    required = {"leftmargin", "rightmargin", "topmargin", "bottommargin",
                "ideogramheight", "ideogramlateralmargin",
                "data1height", "data1inmargin", "data1outmargin",
                "data1min", "data1max",
                "data2height", "data2inmargin", "data2outmargin",
                "data2min", "data2max",
                "dataideogrammin", "dataideogrammax",
                "dataallmin", "dataallmax"}
    assert required.issubset(pp.keys())


def test_invalid_plot_type():
    with pytest.raises(ValueError):
        get_default_plot_params(99)


def test_plot_type_1_specific_values():
    pp = get_default_plot_params(1)
    assert pp["data2height"] == 0
    assert pp["data1height"] == 200
    assert pp["ideogramheight"] == 50
