from numba import jit

@jit(nopython=True, nogil=True)
def remove_peak_above_precursormz(peak_array, precursormz):
    """
    Remove peaks that are above a specified precursor m/z value.

    :param peak_array: A numpy array containing peak data. It's expected to include the m/z (mass-to-charge ratio)
                       values in the first column and corresponding intensities in the second column.

    :param precursormz: A float representing the precursor m/z value. Peaks with m/z values greater than
                        (precursormz + 5 Da) are considered as noise and removed from the peak array.

    :return: The filtered numpy array that only includes peaks with m/z values below the (precursor m/z + 5 Da).
    """
    # check if the input precursormz is a floating-point number
    if isinstance(precursormz, float):
        # filter the array to only include peaks with m/z less than (precursor m/z + 5 Da)
        peak_array = peak_array[peak_array[:, 0] < precursormz + 5.0]  # Filtering criterion
        return peak_array  # Return the filtered array

    return peak_array  # If precursormz is not a float, return the original peak array without any filtering
