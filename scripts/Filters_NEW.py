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

def check_minimum_peak_requiered(peak_dataframe):
    """
    Checks if the given peak_dataframe has the minimum required number of peaks.

    :param peak_dataframe: A dataframe containing peak information.
    :type peak_dataframe: pandas.DataFrame

    :return: The peak_dataframe if it has at least 3 peaks, None otherwise.
    :rtype: pandas.DataFrame or None
    """
    if len(peak_dataframe) < 3:
        return None
    else:
        return peak_dataframe

def normalize_intensity(peak_dataframe):
    """
    Method to normalize the intensity values in a peak dataframe.

    :param peak_dataframe: pandas DataFrame containing peak data with 'intensity' column
    :return: pandas DataFrame with normalized intensity values

    This method uses the formula (value - min) / (max - min) to normalize the intensity values
    in the peak dataframe. The 'intensity' column is modified in-place.
    """
    peak_dataframe['intensity'] = peak_dataframe['intensity'] = peak_dataframe['intensity']/peak_dataframe['intensity'].max()

    return peak_dataframe

def keep_mz_in_range(peak_dataframe):
    """
    Filter the peak dataframe to keep only the rows where the 'mz' column value is within the range of 50 and 2000.

    :param peak_dataframe: A pandas DataFrame containing peaks information.
    :return: A filtered pandas DataFrame with the 'mz' values within the specified range.
    """
    peak_dataframe = peak_dataframe.loc[(peak_dataframe['mz'] >= 50) & (peak_dataframe['mz'] <= 2000)]

    return peak_dataframe

def check_minimum_of_high_peaks_requiered(peak_dataframe, intensity_percent, no_peaks):
    """
    Check if the number of high peaks in a peak dataframe meets the minimum requirement.

    :param peak_dataframe: DataFrame object containing peak data.
    :param intensity_percent: Intensity percentage threshold for considering peaks as high peaks.
    :param no_peaks: Minimum number of high peaks required.

    :return: DataFrame object containing the high peaks if the number of high peaks meets the minimum requirement,
             None otherwise.
    """
    percent_of_max = peak_dataframe['intensity']/peak_dataframe['intensity'].max() * 100
    filtered_df = peak_dataframe[percent_of_max >= intensity_percent]
    if len(filtered_df) < no_peaks:
        return None
    else:
        return filtered_df

def apply_filters(peak_dataframe, precursormz):
    """

    Apply filters on a peak dataframe based on the given parameters.

    :param peak_dataframe: A dataframe containing peak information.
    :param precursormz: The precursor m/z value.
    :return: The filtered peak dataframe.

    """
    peak_dataframe = check_minimum_peak_requiered(peak_dataframe)
    if peak_dataframe is None:
        return None
    else:
        peak_dataframe = remove_peak_above_precursormz(peak_dataframe, precursormz)
        peak_dataframe = reduce_peak_list(peak_dataframe)
        peak_dataframe = normalize_intensity(peak_dataframe)
        peak_dataframe = keep_mz_in_range(peak_dataframe)
        peak_dataframe = check_minimum_of_high_peaks_requiered(peak_dataframe, intensity_percent=5.0, no_peaks=2)
        if peak_dataframe is None:
            return None

        return peak_dataframe