from numba import jit

@jit(nopython=True, nogil=True)
def keep_mz_in_range(peak_array, mz_from, mz_to):
    """
    This function is designed to filter a given peak array based on a specified
    mass-to-charge (m/z) values range. Peaks out of this range will be excluded.

    Arguments:
    :param peak_array: A 2D numpy array where each row represents a peak and each column represents a data field.
                       The first column is assumed to contain the m/z values for the peaks.
    :param mz_from: The lower limit of the m/z range allowed.
    :param mz_to: The upper limit of the m/z range allowed.

    :return: Returns a new numpy array containing only the peaks within the specified m/z range.

    """

    # Define a boolean mask where True values correspond to m/z values within the specified range
    # Here, we directly compare the m/z values in the peak array (assumed to be in the first column)
    # with the specified lower and upper limits.
    # The '&' operator is used to logically 'and' the two condition arrays, so we only get True
    # where both conditions are met
    mz_range = (peak_array[:, 0] >= int(mz_from)) & (peak_array[:, 0] <= int(mz_to))

    # Apply the mask to the peak array to remove peaks outside the specified range, returning this filtered array
    return peak_array[mz_range]
