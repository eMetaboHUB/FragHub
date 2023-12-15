import pandas as pd

def remove_peak_above_precursormz(peak_dataframe, precursormz):
    """
    :param peak_dataframe: The dataframe containing peak information.
    :param precursormz: The precursor m/z value used as a threshold for filtering peaks.
    :return: The filtered dataframe with peaks having m/z values below the provided precursor m/z + 5.0 Da.
    """
    peak_dataframe = peak_dataframe[peak_dataframe['mz'] < float(precursormz) + 5.0]  # Removing peaks with mz > precursormz + 5 Da

    return peak_dataframe

def apply_filters(peak_dataframe, precursormz):
    """
    Apply Filters

    Applies filters to a peak dataframe based on a given precursor m/z.

    Parameters:
    :param peak_dataframe: The dataframe containing peak data.
    :param precursormz: The precursor m/z value to filter the peaks.

    Returns:
    The filtered peak dataframe.

    """
    peak_dataframe = remove_peak_above_precursormz(peak_dataframe, precursormz)

    return peak_dataframe