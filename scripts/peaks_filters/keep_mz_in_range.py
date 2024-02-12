
def keep_mz_in_range(peak_array, mz_from, mz_to):
    """
    :param peak_array: A numpy array containing peak data. The array must have two columns with the first column for 'mz'.
    :param mz_from: The lower bound of the m/z range. Peaks with m/z values less than this threshold will be excluded from the filtered array.
    :param mz_to: The upper bound of the m/z range. Peaks with m/z values greater than this threshold will be excluded from the filtered array.

    :return: A new numpy array containing only the peaks within the specified m/z range.
    """
    mz_range = (peak_array[:,0] >= int(mz_from)) & (peak_array[:,0] <= int(mz_to))
    return peak_array[mz_range]
