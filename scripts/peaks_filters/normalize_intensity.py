import numpy as np

def normalize_intensity(peak_array):
    """
    Normalize (rescale) the intensity values of a numpy array.

    :param peak_array: A numpy array containing peak data.
                       The array is assumed to have two columns,
                       with the second column containing the intensities.
    :type peak_array: numpy.ndarray
    :return: The input numpy array with normalized intensity values.
    :rtype: numpy.ndarray
    """
    if peak_array[:, 1].max() != 0:
        peak_array[:, 1] = peak_array[:, 1] / peak_array[:, 1].max()

        return peak_array

    return np.array([])
