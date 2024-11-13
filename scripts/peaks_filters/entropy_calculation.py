import numpy as np

def entropy_calculation(peak_list):
    """
    Calculate the Shannon entropy of the peak list.

    :param peak_list: List or numpy array containing peak intensities
    :return: Calculated entropy value
    """
    # Convert peak list to a numpy array if it's not already
    if not isinstance(peak_list, np.ndarray):
        peak_list = np.array(peak_list)

    # Normalize the peak list to get probabilities
    total_intensity = np.sum(peak_list)
    if total_intensity == 0:
        return 0
    probabilities = peak_list / total_intensity

    # Filter out zero probabilities to avoid log(0)
    probabilities = probabilities[probabilities > 0]

    # Calculate the Shannon entropy
    entropy = -np.sum(probabilities * np.log2(probabilities))

    return entropy