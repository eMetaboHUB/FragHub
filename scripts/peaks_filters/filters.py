from .check_minimum_of_high_peaks_requiered import *
from .remove_peak_above_precursormz import *
from .check_minimum_peak_requiered import *
from .normalize_intensity import *
from .keep_mz_in_range import *
from .reduce_peak_list import *
import numpy as np

def apply_filters(peak_array, precursormz, parameters_dict):
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
        if peak_array.size == 0:
            return np.array([])

    if parameters_dict['remove_peak_above_precursormz'] == 1.0:
        peak_array = remove_peak_above_precursormz(peak_array, precursormz)
        if peak_array.size == 0:
            return np.array([])

    if parameters_dict['reduce_peak_list'] == 1.0:
        peak_array = reduce_peak_list(peak_array, max_peaks)

    if parameters_dict['normalize_intensity'] == 1.0:
        peak_array = normalize_intensity(peak_array)
        if peak_array.size == 0:
            return np.array([])

    if parameters_dict['keep_mz_in_range'] == 1.0:
        peak_array = keep_mz_in_range(peak_array, mz_from, mz_to)

    if parameters_dict['check_minimum_of_high_peaks_requiered'] == 1.0:
        peak_array = check_minimum_of_high_peaks_requiered(peak_array, intensity_percent, no_peaks)

    if peak_array.size == 0:
        return np.array([])

    return peak_array
