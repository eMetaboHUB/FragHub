from set_parameters import parameters_dict
import pandas as pd
import numpy as np

def remove_peak_above_precursormz(peak_array, precursormz):
    """
    :param peak_array: numpy array containing peak data
    :param precursormz: float representing the precursor m/z value
    :return: filtered numpy array with peaks below the specified precursor m/z value + 5 Da
    """
    if not isinstance(precursormz, float):
        try:
            precursormz = float(precursormz)
            peak_array = peak_array[peak_array[:,0] < precursormz + 5.0]  # Removing peaks with mz > precursormz + 5 Da
            return peak_array
        except:
            return np.empty((0,2))

    return peak_array

def reduce_peak_list(peak_array, max_peaks):
    """
    Reduce the peak list to a specified number of maximum peaks.

    :param peak_array: The numpy array of peak data.
    :param max_peaks: The maximum number of peaks to retain.
    :return: The reduced peak numpy array.
    """
    if len(peak_array) > int(max_peaks):
        # Sort the array by intensity in descending order and keep the top rows
        peak_array = peak_array[peak_array[:,1].argsort()[::-1][:int(max_peaks)]]
        # Then, re-sort by 'mz'
        peak_array = peak_array[peak_array[:,0].argsort()]

    return peak_array

def check_minimum_peak_requiered(peak_array, n_peaks):
    """
    :param peak_array: A numpy array representing the peaks
    :param n_peaks: An integer representing the minimum number of peaks required
    :return: A numpy array representing the peaks if the number of peaks is not less than n_peaks.
             Otherwise, an empty array is returned.
    """
    if len(peak_array) < int(n_peaks):
        return np.empty((0,2))
    else:
        return peak_array

def normalize_intensity(peak_array):
    """
    Normalize the intensity values of a numpy array.

    :param peak_array: A numpy array containing peak data.
                       The array is assumed to have two columns,
                       with the second column containing the intensities.
    :type peak_array: numpy.ndarray
    :return: The input numpy array with normalized intensity values.
    :rtype: numpy.ndarray
    """
    peak_array[:, 1] = peak_array[:, 1] / peak_array[:, 1].max()

    return peak_array

def keep_mz_in_range(peak_array, mz_from, mz_to):
    """
    :param peak_array: A numpy array containing peak data. The array must have two columns with the first column for 'mz'.
    :param mz_from: The lower bound of the m/z range. Peaks with m/z values less than this threshold will be excluded from the filtered array.
    :param mz_to: The upper bound of the m/z range. Peaks with m/z values greater than this threshold will be excluded from the filtered array.

    :return: A new numpy array containing only the peaks within the specified m/z range.
    """
    mz_range = (peak_array[:,0] >= int(mz_from)) & (peak_array[:,0] <= int(mz_to))
    return peak_array[mz_range]

def check_minimum_of_high_peaks_requiered(peak_array, intensity_percent, no_peaks): # NOTE: potentiel probleme avec 3659e269-2355-485c-bfd8-cacc2a488a3e, retourn []
    """
    :param peak_array: A numpy array containing peak data. The array must have two columns with the first column for 'intensity'.
    :param intensity_percent: The minimum percentage of the maximum intensity required for a peak to be considered high.
    :param no_peaks: The minimum number of high peaks required.

    :return: If the number of high peaks in peak_array is less than no_peaks, an empty numpy array is returned. Otherwise, a filtered array containing the high peaks is returned.
    """
    if peak_array.size == 0:
        return peak_array

    percent_of_max = (peak_array[:,1] / peak_array[:,1].max()) * 100
    filtered_array = peak_array[percent_of_max >= intensity_percent]

    if len(filtered_array) < int(no_peaks):
        return np.empty((0,2))
    else:
        return filtered_array

def apply_filters(peak_array, precursormz):
    """
    :param peak_array: the input numpy array containing peak information
    :param precursormz: the precursor m/z value
    :return: a filtered numpy array containing peak information

    """
    n_peaks = parameters_dict['check_minimum_peak_requiered_n_peaks']
    max_peaks = parameters_dict['reduce_peak_list_max_peaks']
    mz_from = parameters_dict['keep_mz_in_range_from_mz']
    mz_to = parameters_dict['keep_mz_in_range_to_mz']
    intensity_percent = parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent']
    no_peaks = parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks']

    if parameters_dict['check_minimum_peak_requiered'] == 1.0:
        peak_array = check_minimum_peak_requiered(peak_array, n_peaks)
        if len(peak_array) == 0:
            return np.array([])

    if parameters_dict['remove_peak_above_precursormz'] == 1.0:
        peak_array = remove_peak_above_precursormz(peak_array, precursormz)
        if len(peak_array) == 0:
            return np.array([])

    if parameters_dict['reduce_peak_list'] == 1.0:
        peak_array = reduce_peak_list(peak_array, max_peaks)

    if parameters_dict['normalize_intensity'] == 1.0:
        peak_array = normalize_intensity(peak_array)

    if parameters_dict['keep_mz_in_range'] == 1.0:
        peak_array = keep_mz_in_range(peak_array, mz_from, mz_to)

    if parameters_dict['check_minimum_of_high_peaks_requiered'] == 1.0:
        peak_array = check_minimum_of_high_peaks_requiered(peak_array, intensity_percent, no_peaks)

    if len(peak_array) == 0:
        return np.array([])

    return peak_array