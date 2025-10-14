from numba import jit
import numpy as np


@jit(nopython=True, nogil=True)
def normalize_intensity(peak_array):
    """
    Normalize (rescale) the intensity values of a numpy array.

    This function checks if the maximum value of intensity (2nd column of the array) is non-zero,
    if it is then the function divides all the intensity values by the maximum intensity value
    to normalize them. If the maximum intensity is zero, the function returns an empty numpy array.

    :param peak_array: A numpy array containing peak data.
                       The array is assumed to have two columns,
                       with the second column containing the intensities.
    :type peak_array: numpy.ndarray
    :return: The input numpy array with normalized intensity values.
             If the max intensity is zero, an empty array is returned.
    :rtype: numpy.ndarray
    """
    # Check if the maximum intensity value is not zero
    if peak_array[:, 1].max() != 0:
        # If so, normalize all intensity values by the maximum value
        peak_array[:, 1] = peak_array[:, 1] / peak_array[:, 1].max()
        # Return the array with normalized intensity
        return peak_array

    # If maximum intensity is zero, return an empty array
    return np.empty((0, 2), dtype=np.float64)
