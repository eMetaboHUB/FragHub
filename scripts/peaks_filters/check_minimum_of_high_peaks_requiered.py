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
    # The guard clause check if the input peak array is empty
    if peak_array.size == 0:
        return peak_array

    # Calculate the percentage of maximum for each peak intensity in the array
    percent_of_max = (peak_array[:, 1] / peak_array[:, 1].max()) * 100

    # Filter the peak array based on the intensity percentage requirement
    filtered_array = peak_array[percent_of_max >= intensity_percent]

    # If the number of peaks in the filtered array is less than the required number of peaks,
    # return an empty array (0,2)
    if len(filtered_array) < int(no_peaks):
        return np.empty((0, 2))
    else:  # Otherwise, return the original peak array
        return peak_array
