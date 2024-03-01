
def remove_peak_above_precursormz(peak_array, precursormz):
    """
    :param peak_array: numpy array containing peak data
    :param precursormz: float representing the precursor m/z value
    :return: filtered numpy array with peaks below the specified precursor m/z value + 5 Da
    """
    if isinstance(precursormz, float):
        peak_array = peak_array[peak_array[:,0] < precursormz + 5.0]  # Removing peaks with mz > precursormz + 5 Da
        return peak_array

    return peak_array
