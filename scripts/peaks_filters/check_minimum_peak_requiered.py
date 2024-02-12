import numpy as np

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
