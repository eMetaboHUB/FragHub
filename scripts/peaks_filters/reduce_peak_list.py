from numba import jit

@jit(nopython=True, nogil=True)
def reduce_peak_list(peak_array, max_peaks):
    """
    Reduce the peak list to a specified number of maximum peaks.

    :param peak_array: The numpy array of peak data.
    :param max_peaks: The maximum number of peaks to retain.
    :return: The reduced peak numpy array.
    """
    # Check if the length of the peak array is greater than the maximum number of peaks allowed
    if len(peak_array) > int(max_peaks):
        # If it is, sort the array by the second column (intensity values) in descending order (highest to lowest)
        # The sorted function returns indices that would sort the array
        # We then use these indices to get the top 'max_peaks' rows from the array. This is done by slicing the
        # indices array from the end towards the start up to 'max_peaks' and reverse it (to get descending order)
        peak_array = peak_array[peak_array[:, 1].argsort()[::-1][:int(max_peaks)]]
        # After the intensity-based sorting and slicing, we then sort these selected peaks by their 'mz' values
        # (first column in the array). We again use the argsort function, but this time without reversing the array
        # as we want 'mz' to be sorted in ascending order
        peak_array = peak_array[peak_array[:, 0].argsort()]
    # After applying the objective conditions, we return the final reduced array of peaks
    return peak_array
