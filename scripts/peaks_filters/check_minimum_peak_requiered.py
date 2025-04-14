import scripts.deletion_report
import numpy as np

def check_minimum_peak_requiered(spectrum, peak_array, n_peaks):
    """
    This function checks if the number of peaks in an array is not less than a certain amount.

    :param peak_array: A numpy array representing the peaks. This array holds all the peak
    information that will be checked against our minimum number of peaks (n_peaks).

    :param n_peaks: An integer representing the minimum number of peaks required. If the number
    of peaks in 'peak_array' is less than this value, an empty array will be returned.


    :return: If the number of peaks in 'peak_array' is not less than 'n_peaks', 'peak_array' is
    returned as is. But, if 'peak_array' has fewer peaks than 'n_peaks', an empty numpy array is
    returned. This is represented as a (0,2) shape numpy array, indicating no data in two dimensions.
    """
    if len(peak_array) < int(n_peaks):  # If the amount of peaks in 'peak_array' is less than 'n_peaks'...
        spectrum['DELETION_REASON'] = "spectrum deleted because its number of peaks is below the threshold chosen by the user"
        scripts.deletion_report.deleted_spectrum_list.append(spectrum)
        scripts.deletion_report.minimum_peaks_not_requiered += 1
        return np.empty((0, 2))  # return an empty numpy array
    else:
        return peak_array  # Otherwise, return 'peak_array' as is
