from pyteomics import achrom
import matplotlib.pyplot as plt

def predict_lc_retention_times(peptides):
    """
    Map each peptide sequence to its predicted retention time.

    Parameters
    ----------
    peptides : list of str
        A list containing the sequences of the peptides.
    Returns
    -------
    dict
        With peptide sequences as keys and predicted retention times (rounded to 2 decimals) as values.
    """

    lc_retention_times = {}
    for peptide in peptides:
        peptide_RT = achrom.calculate_RT(peptide, achrom.RCs_guo_ph7_0)
        lc_retention_times[peptide] = round(peptide_RT, 2)

    return lc_retention_times


def plot_retention_time(retention_times, resolution=30):
    """
    Create a histogram of the retention times.

    Parameters
    ----------
    retention_times : list|dict
        A list or dictionary with the retention times.
    resolution : int, optional
        The number of histogram bins.
    Returns
    -------
    None
    """

    if type(retention_times) == list:
        plt.hist(retention_times, bins=resolution)
    elif type(retention_times) == dict:
        plt.hist(retention_times.values(), bins=resolution)

    plt.xlabel("retention time (minutes)")
    plt.ylabel("number of peptides")
    plt.show()


def select_retention_time_window(peptide_rt_map, lower_ret_time, upper_ret_time):
    """
    Accept a dictionary of peptide sequences and their retention times and return only the peptides with retention times inside
    the given lower and upper bounds. Returns a dictionary with the peptides and their retention times.

    Parameters
    ----------
    peptide_rt_map : dict
        With peptide sequences as keys and retention times as values.
    lower_ret_time : float
        Lower bound for the retention time.
    upper_ret_time: float
        Upper bound for the retention time.
    Returns
    -------
    dict
        With filtered peptide sequences as keys and retention times as values.
    """

    filtered_rt_map = {}
    for peptide, rt in peptide_rt_map.items():
        if lower_ret_time <= rt <= upper_ret_time:
            filtered_rt_map[peptide] = rt

    return filtered_rt_map