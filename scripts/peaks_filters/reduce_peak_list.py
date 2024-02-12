
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
