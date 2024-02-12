import numpy as np

def check_minimum_of_high_peaks_requiered(peak_array, intensity_percent, no_peaks):
    """
    :param peak_array: An array containing peak values and intensities.
    :param intensity_percent: The minimum percentage of maximum intensity required for a peak to be considered.
    :param no_peaks: The minimum number of peaks required.

    :return: If the peak_array is empty, it returns peak_array.
             If the number of peaks in filtered_array is less than no_peaks, it returns an empty array (0,2).
             Otherwise, it returns peak_array.

    """
    if peak_array.size == 0:
        return peak_array

    percent_of_max = (peak_array[:,1] / peak_array[:,1].max()) * 100
    filtered_array = peak_array[percent_of_max >= intensity_percent]

    if len(filtered_array) < int(no_peaks):
        return np.empty((0,2))
    else:
        return peak_array
