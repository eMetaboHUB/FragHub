from .check_minimum_of_high_peaks_requiered import *
from .remove_peak_above_precursormz import *
from .check_minimum_peak_requiered import *
from .normalize_intensity import *
from .keep_mz_in_range import *
from .reduce_peak_list import *
import numpy as np

def apply_filters(peak_array, precursormz, parameters_dict, peaks_filters_not_valid):
    """
    Function to apply various filters on a given peak_array according to the provided parameters.

    :param peak_array: the input numpy array containing peak information
    :param precursormz: the precursor m/z value
    :param parameters_dict: a dictionary containing various filtering parameters
    :return: a filtered numpy array containing peak information
    """
    # retrieve filtering parameters from the dictionary
    n_peaks = parameters_dict['check_minimum_peak_requiered_n_peaks']
    max_peaks = parameters_dict['reduce_peak_list_max_peaks']
    mz_from = parameters_dict['keep_mz_in_range_from_mz']
    mz_to = parameters_dict['keep_mz_in_range_to_mz']
    intensity_percent = parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent']
    no_peaks = parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks']

    # apply filters in a sequence
    if parameters_dict['check_minimum_peak_requiered'] == 1.0:
        # filter out peaks below a minimum threshold
        peak_array = check_minimum_peak_requiered(peak_array, n_peaks)
        # if no peaks pass this filter, return an empty array
        if peak_array.size == 0:
            peaks_filters_not_valid = 1
            return np.array([]), peaks_filters_not_valid

    if parameters_dict['remove_peak_above_precursormz'] == 1.0:
        # remove peaks that are above a specified limit
        peak_array = remove_peak_above_precursormz(peak_array, precursormz)
        # if no peaks pass this filter, return an empty array
        if peak_array.size == 0:
            peaks_filters_not_valid = 1
            return np.array([]), peaks_filters_not_valid

    if parameters_dict['reduce_peak_list'] == 1.0:
        # limit total number of peaks to be considered
        peak_array = reduce_peak_list(peak_array, max_peaks)

    if parameters_dict['normalize_intensity'] == 1.0:
        # normalize the intensity of peaks
        peak_array = normalize_intensity(peak_array)
        # if no peaks pass this filter, return an empty array
        if peak_array.size == 0:
            peaks_filters_not_valid = 1
            return np.array([]), peaks_filters_not_valid

    if parameters_dict['keep_mz_in_range'] == 1.0:
        # keep peaks within a certain mz range
        peak_array = keep_mz_in_range(peak_array, mz_from, mz_to)

    if parameters_dict['check_minimum_of_high_peaks_requiered'] == 1.0:
        # filter out peaks below a minimum intensity percent
        peak_array = check_minimum_of_high_peaks_requiered(peak_array, intensity_percent, no_peaks)

    # if no peaks pass the filters, return an empty array
    if peak_array.size == 0:
        peaks_filters_not_valid = 1
        return np.array([]), peaks_filters_not_valid

    # return the filtered peak array
    return peak_array, peaks_filters_not_valid
