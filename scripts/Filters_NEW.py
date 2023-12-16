from set_parameters import parameters_dict
import pandas as pd

def remove_peak_above_precursormz(peak_dataframe, precursormz):
    """
    :param peak_dataframe: pandas DataFrame containing peak data
    :param precursormz: float representing the precursor m/z value
    :return: filtered pandas DataFrame with peaks below the specified precursor m/z value + 5 Da

    """
    peak_dataframe = peak_dataframe[peak_dataframe['mz'] < float(precursormz) + 5.0]  # Removing peaks with mz > precursormz + 5 Da

    return peak_dataframe

def reduce_peak_list(peak_dataframe,max_peaks):
    """
    Reduce the peak list to a specified number of maximum peaks.

    :param peak_dataframe: The dataframe of peak data.
    :param max_peaks: The maximum number of peaks to retain.
    :return: The reduced peak dataframe.
    """
    if len(peak_dataframe) > int(max_peaks):
        # Triez le DataFrame par intensité en ordre décroissant et gardez les 500 premières lignes
        peak_dataframe = peak_dataframe.sort_values('intensity', ascending=False).head(int(max_peaks))
        # Re-classez ensuite par 'mz' si vous le souhaitez
        peak_dataframe = peak_dataframe.sort_values('mz')

    return peak_dataframe

def check_minimum_peak_requiered(peak_dataframe,n_peaks):
    """
    :param peak_dataframe: A pandas DataFrame representing the peaks
    :param n_peaks: An integer representing the minimum number of peaks required
    :return: A pandas DataFrame representing the peaks if the number of peaks is not less than n_peaks.
             Otherwise, an empty DataFrame is returned.
    """
    if len(peak_dataframe) < int(n_peaks):
        return pd.DataFrame()
    else:
        return peak_dataframe

def normalize_intensity(peak_dataframe):
    """
    Normalize the intensity values of a DataFrame.

    :param peak_dataframe: A pandas DataFrame containing peak data.
    :type peak_dataframe: pandas.DataFrame
    :return: The input DataFrame with normalized intensity values.
    :rtype: pandas.DataFrame
    """
    peak_dataframe.loc[:, 'intensity'] = peak_dataframe['intensity']/peak_dataframe['intensity'].max()

    return peak_dataframe

def keep_mz_in_range(peak_dataframe,mz_from,mz_to):
    """

    :param peak_dataframe: A pandas DataFrame containing peak data. Each row represents a peak, and the DataFrame must have a column named 'mz' containing the mass-to-charge ratio (m/z)
    * values of the peaks.

    :param mz_from: The lower bound of the m/z range. Peaks with m/z values less than this threshold will be excluded from the filtered DataFrame.

    :param mz_to: The upper bound of the m/z range. Peaks with m/z values greater than this threshold will be excluded from the filtered DataFrame.

    :return: A new pandas DataFrame containing only the peaks within the specified m/z range.

    """
    peak_dataframe = peak_dataframe.loc[(peak_dataframe['mz'] >= int(mz_from)) & (peak_dataframe['mz'] <= int(mz_to))]

    return peak_dataframe

def check_minimum_of_high_peaks_requiered(peak_dataframe, intensity_percent, no_peaks):
    """
    :param peak_dataframe: A pandas DataFrame containing peak data. The DataFrame must have columns named 'intensity' and 'intensity_percent'.
    :param intensity_percent: The minimum percentage threshold of the maximum intensity required for a peak to be considered high.
    :param no_peaks: The minimum number of high peaks required.

    :return: If the number of high peaks in peak_dataframe is less than no_peaks, an empty pandas DataFrame is returned. Otherwise, a filtered DataFrame containing the high peaks is returned
    *.
    """
    percent_of_max = peak_dataframe['intensity']/peak_dataframe['intensity'].max() * 100
    filtered_df = peak_dataframe[percent_of_max >= intensity_percent]
    if len(filtered_df) < int(no_peaks):
        return pd.DataFrame()
    else:
        return filtered_df

def apply_filters(peak_dataframe, precursormz):
    """
    :param peak_dataframe: the input dataframe containing peak information
    :param precursormz: the precursor m/z value
    :return: a filtered dataframe containing peak information

    """
    n_peaks = parameters_dict['check_minimum_peak_requiered_n_peaks']
    max_peaks = parameters_dict['reduce_peak_list_max_peaks']
    mz_from = parameters_dict['keep_mz_in_range_from_mz']
    mz_to = parameters_dict['keep_mz_in_range_to_mz']
    intensity_percent = parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent']
    no_peaks = parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks']


    peak_dataframe = check_minimum_peak_requiered(peak_dataframe,n_peaks)
    if peak_dataframe.empty:
        return pd.DataFrame()
    else:
        if parameters_dict['remove_peak_above_precursormz'] == 1.0:
            peak_dataframe = remove_peak_above_precursormz(peak_dataframe, precursormz)
        if parameters_dict['reduce_peak_list'] == 1.0:
            peak_dataframe = reduce_peak_list(peak_dataframe,max_peaks)
        if parameters_dict['normalize_intensity'] == 1.0:
            peak_dataframe = normalize_intensity(peak_dataframe)
        if parameters_dict['keep_mz_in_range'] == 1.0:
            peak_dataframe = keep_mz_in_range(peak_dataframe,mz_from,mz_to)
        if parameters_dict['check_minimum_of_high_peaks_requiered'] == 1.0:
            peak_dataframe = check_minimum_of_high_peaks_requiered(peak_dataframe, intensity_percent, no_peaks)
        if peak_dataframe.empty:
            return pd.DataFrame()

        return peak_dataframe