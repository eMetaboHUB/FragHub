import pandas as pd

def remove_peak_above_precursormz(peak_dataframe, precursormz):
    """
    :param peak_dataframe: The dataframe containing peak information.
    :param precursormz: The precursor m/z value used as a threshold for filtering peaks.
    :return: The filtered dataframe with peaks having m/z values below the provided precursor m/z + 5.0 Da.
    """
    peak_dataframe = peak_dataframe[peak_dataframe['mz'] < float(precursormz) + 5.0]  # Removing peaks with mz > precursormz + 5 Da

    return peak_dataframe

def reduce_peak_list(peak_dataframe):
    """
    :param peak_dataframe: A DataFrame containing peak data.
    :return: A reduced DataFrame containing the top 500 peaks based on intensity, sorted by mz.

    The reduce_peak_list method takes a peak_dataframe as input and returns a reduced DataFrame containing the top 500 peaks based on intensity. The method first sorts the peak_dataframe
    * in descending order based on intensity and retains the top 500 rows. It then sorts the resulting DataFrame by mz.
    """
    if len(peak_dataframe) > 500:
        # Triez le DataFrame par intensité en ordre décroissant et gardez les 500 premières lignes
        peak_dataframe = peak_dataframe.sort_values('intensity', ascending=False).head(500)
        # Re-classez ensuite par 'mz' si vous le souhaitez
        peak_dataframe = peak_dataframe.sort_values('mz')

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
    peak_dataframe = reduce_peak_list(peak_dataframe)

    return peak_dataframe