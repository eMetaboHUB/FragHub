import pandas as pd

def remove_peak_above_precursormz(peak_dataframe, precursormz):
    """
    Removes peaks from `peak_dataframe` with mz values greater than `precursormz` + 5.

    :param peak_dataframe: A DataFrame containing peak data.
    :param precursormz: The precursor m/z value.
    :return: A DataFrame with peaks removed.
    """
    peak_dataframe = peak_dataframe[peak_dataframe['mz'] < float(precursormz) + 5.0]  # Removing peaks with mz > precursormz + 5 Da

    return peak_dataframe

def reduce_peak_list(peak_dataframe):
    """
    Reduce the peak list to a maximum of 500 peaks.

    :param peak_dataframe: The dataframe containing the peak data.
    :return: The reduced peak dataframe.
    """
    if len(peak_dataframe) > 500:
        # Triez le DataFrame par intensité en ordre décroissant et gardez les 500 premières lignes
        peak_dataframe = peak_dataframe.sort_values('intensity', ascending=False).head(500)
        # Re-classez ensuite par 'mz' si vous le souhaitez
        peak_dataframe = peak_dataframe.sort_values('mz')

    return peak_dataframe

def check_minimum_peak_requiered(peak_dataframe):
    """
    Check if the given peak dataframe meets the minimum requirement.

    :param peak_dataframe: A pandas dataframe representing the peak data.
    :return: A pandas dataframe. If the length of the peak dataframe is less than 3, an empty dataframe is returned. Otherwise, the peak dataframe is returned as is.
    """
    if len(peak_dataframe) < 3:
        return pd.DataFrame()
    else:
        return peak_dataframe

def normalize_intensity(peak_dataframe):
    """
    Normalize the intensity values in the given peak_dataframe.

    :param peak_dataframe: The dataframe containing peak information.
    :type peak_dataframe: pandas.DataFrame
    :return: The normalized peak_dataframe.
    :rtype: pandas.DataFrame
    """
    peak_dataframe.loc[:, 'intensity'] = peak_dataframe['intensity']/peak_dataframe['intensity'].max()

    return peak_dataframe

def keep_mz_in_range(peak_dataframe):
    """
        Keep M/Z in Range

        Keeps the M/Z values within a specified range in the given peak dataframe.

        :param peak_dataframe: A pandas DataFrame containing peak data.
        :return: A filtered pandas DataFrame with M/Z values within the specified range.
    """
    peak_dataframe = peak_dataframe.loc[(peak_dataframe['mz'] >= 50) & (peak_dataframe['mz'] <= 2000)]

    return peak_dataframe

def check_minimum_of_high_peaks_requiered(peak_dataframe, intensity_percent, no_peaks):
    """
    Checks if the minimum number of high peaks is met.

    :param peak_dataframe: DataFrame containing peak data.
    :param intensity_percent: Minimum intensity percentage for a peak to be considered high.
    :param no_peaks: Minimum number of high peaks required.
    :return: DataFrame containing the filtered high peaks if the minimum is met, otherwise an empty DataFrame.
    """
    percent_of_max = peak_dataframe['intensity']/peak_dataframe['intensity'].max() * 100
    filtered_df = peak_dataframe[percent_of_max >= intensity_percent]
    if len(filtered_df) < no_peaks:
        return pd.DataFrame()
    else:
        return filtered_df

def apply_filters(peak_dataframe, precursormz):
    """
    Apply various filters to a peak dataframe.

    :param peak_dataframe: A dataframe containing peak information.
    :param precursormz: The precursor m/z value.
    :return: A modified peak dataframe after applying filters.
    """
    peak_dataframe = check_minimum_peak_requiered(peak_dataframe)
    if peak_dataframe.empty:
        return pd.DataFrame()
    else:
        peak_dataframe = remove_peak_above_precursormz(peak_dataframe, precursormz)
        peak_dataframe = reduce_peak_list(peak_dataframe)
        peak_dataframe = normalize_intensity(peak_dataframe)
        peak_dataframe = keep_mz_in_range(peak_dataframe)
        peak_dataframe = check_minimum_of_high_peaks_requiered(peak_dataframe, intensity_percent=5.0, no_peaks=2)
        if peak_dataframe.empty:
            return pd.DataFrame()

        return peak_dataframe