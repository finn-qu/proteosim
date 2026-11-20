from pyteomics import achrom
import matplotlib.pyplot as plt

from proteosim.liquid_chromatography import (
    predict_lc_retention_times,
    plot_retention_time,
    select_retention_time_window
)

def test_predict_lc_retention_times():
    peptides = ["KPRWWWK", "VYTKRQ"]
    expected = {
        peptides[0]: round(achrom.calculate_RT(peptides[0], achrom.RCs_guo_ph7_0), 2),
        peptides[1]: round(achrom.calculate_RT(peptides[1], achrom.RCs_guo_ph7_0), 2)
        }

    actual = predict_lc_retention_times(peptides)
    assert actual == expected


def test_select_retention_time_window():
    peptide_rt_map = {
        "A": 1.0,
        "K": 10.0,
        "V": 100.0
        }
    selected1 = select_retention_time_window(peptide_rt_map, lower_ret_time=20, upper_ret_time=50)      # should return an empty dictionary
    selected2 = select_retention_time_window(peptide_rt_map, lower_ret_time=1, upper_ret_time=100)      # should return 1, 10 and 100
    selected3 = select_retention_time_window(peptide_rt_map, lower_ret_time=0.5, upper_ret_time=100.5)  # should return 1, 10 and 100
    selected4 = select_retention_time_window(peptide_rt_map, lower_ret_time=2, upper_ret_time=99)       # should return 10

    assert selected1 == {
        }
    assert selected2 == {
        "A": 1.0,
        "K": 10.0,
        "V": 100.0
        }
    assert selected3 == {
        "A": 1.0,
        "K": 10.0,
        "V": 100.0
        }
    assert selected4 == {
        "K": 10.0
        }
