from scripts.GUI.utils.global_vars import parameters_dict
from scripts.peaks_filters.entropy_calculation import *
from scripts.calculate_maximized_chunk_size import *
from scripts.normalizer.values_normalizer import *
from scripts.peaks_filters.filters import *
import concurrent.futures
import scripts.deletion_report
import scripts.globals_vars
import numpy as np
import re

np.set_printoptions(suppress=True)

def peak_list_cleaning(spectrum, peak_list, precursormz):
    """
    This function serves the purpose of converting a list of peak tuples (mz and intensity values) into a numpy array.

    :param peak_list: A list of tuples. Each tuple represents a peak and it contains mz value and intensity of the peak.
    :param precursormz: The precursor mz value for the peaks.

    :return: The function returns a numpy array that represents the peaks.
             The peaks in the returned numpy array are sorted by their mz values.
             The peaks are also filtered based on the precursor mz value.
    """

    # Convert the list of tuples (peak_list) to a numpy array.
    # Each tuple in the list contains two floating point values - mz value and intensity of a peak.
    peak_list = np.array(peak_list, dtype=float)

    # Sort the numpy array of peaks based on the mz values.
    # The mz values are the first value in each tuple, hence located in the first column of the array.
    # Method np.argsort() returns indices that sort the mz values.
    # By indexing the original array with these indices, we get sorted array.
    peak_list = peak_list[peak_list[:, 0].argsort()]

    # Apply appropriate set of filters to the sorted numpy array of peaks.
    # The filters use the precursor mz value and pre-defined parameters.
    peak_list = apply_filters(spectrum, peak_list, precursormz, parameters_dict)

    # Return the sorted and preset-filtered numpy array of peaks.
    return peak_list

def peak_list_to_str(peak_list_np):
    """
    This function converts a peak list numpy array into a string format. Each peak is represented as a space-separated
    string of floating point values, and each peak is separated by a newline.

    :param peak_list_np: The peak list represented as a numpy array. The numpy array is expected to contain floating-point values.
    :return: A string representation of the peak list where each row represents a peak and is formatted as a space-separated
    string of floating-point values.
    """
    # Convert the peak list numpy array to a list of lists. Each sub-list contains floating point values representing a peak.
    # The numpy round function is used to limit the precision of the floating-point values to 8 decimal places.
    peak_list_np = peak_list_np.round(8).tolist()

    # Convert the list of lists into a string format. Each sub-list (representing a peak) is converted into a space-separated
    # string of floats. Sub-lists are then joined together with newline characters to form the final string.
    #
    # The join function is used twice:
    # 1) The inner join function is used to convert each sub-list into a space-separated string of floats.
    #    The string formatter is used to ensure each float is represented to 8 decimal places.
    # 2) The outer join function is used to join together all the space-separated strings (representing peaks) with newline characters.
    peak_list_np = "\n".join(" ".join(f"{val:.8f}" for val in sublist) for sublist in peak_list_np)

    # Return the final string representation of the peak list.
    return peak_list_np

def spectrum_cleaning(spectrum):
    """
    This function cleans the input 'spectrum'.

    It starts by checking if the 'PEAKS_LIST' attribute of the spectrum exists.
    If it does, it continues with normalizing all key-values in the spectrum.
    After normalization, the function checks if the 'PRECURSORMZ' exists in the spectrum and a regular expression search for
    the 'float_check_pattern' in the 'PRECURSORMZ' key returns a match.
    If both conditions are true, it modifies 'PRECURSORMZ' to match the regex group, converts it to float and replaces any commas with periods.
    If the float value is positive, it proceeds to convert the peak list to a numpy array and then to string,
    subsequently updating the spectrum's 'PEAKS_LIST' and 'NUM PEAKS' attributes.

    :param spectrum: dictionary containing spectrum information
    :return: cleaned spectrum dictionary if it passes all checks, otherwise None
    """
    peak_list = spectrum["PEAKS_LIST"]
    # If peak_list is not present in the spectrum dictionary, it returns None
    if not peak_list:
        spectrum['DELETION_REASON'] = "spectrum deleted because peaks list is empty"
        scripts.deletion_report.deleted_spectrum_list.append(spectrum)
        scripts.deletion_report.no_peaks_list += 1
        return None
    spectrum = normalize_values(spectrum)
    # If normalization of spectrum fails, it returns None
    if not spectrum:
        return None
    # Checks if "PRECURSORMZ" exists in the spectrum
    if "PRECURSORMZ" in spectrum and ("_GC" not in spectrum["FILENAME"] and not re.search(r"\bGC\b", spectrum["INSTRUMENTTYPE"])):
        if re.search(scripts.globals_vars.float_check_pattern, str(spectrum["PRECURSORMZ"])):
            # 'PRECURSORMZ' modification with match from a regular expression search for 'float_check_pattern'
            spectrum["PRECURSORMZ"] = re.search(scripts.globals_vars.float_check_pattern, str(spectrum["PRECURSORMZ"])).group(1)
            float_precursor_mz = float(spectrum["PRECURSORMZ"].replace(",", "."))
            # Float value of 'PRECURSORMZ' needs to be greater than 0
            if float_precursor_mz <= 0.0:
                spectrum['DELETION_REASON'] = "spectrum deleted because precursor mz is less than or equal to zero."
                scripts.deletion_report.deleted_spectrum_list.append(spectrum)
                scripts.deletion_report.no_precursor_mz += 1
                return None
            # Converts peak list to a numpy array
            peak_list_np = peak_list_cleaning(spectrum, peak_list, float_precursor_mz)

            if peak_list_np.size == 0:
                return None

            intensities = peak_list_np[:, 1]
            spectrum["ENTROPY"] = str(entropy_calculation(intensities))
            if parameters_dict["remove_spectrum_under_entropy_score"] == 1.0:
                if re.search(scripts.globals_vars.float_check_pattern, str(spectrum["ENTROPY"])):
                    if float(spectrum["ENTROPY"]) < parameters_dict["remove_spectrum_under_entropy_score_value"]:
                        spectrum['DELETION_REASON'] = "spectrum deleted because it's entropy score is lower than the threshold selected by the user."
                        scripts.deletion_report.deleted_spectrum_list.append(spectrum)
                        scripts.deletion_report.low_entropy_score += 1
                        return None
            # If numpy array is empty, it returns none
            if peak_list_np.size == 0:
                return None
            spectrum["NUM PEAKS"] = str(peak_list_np.shape[0])
            # Convert numpy array back to string and update 'PEAKS_LIST' in spectrum
            peak_list_np = peak_list_to_str(peak_list_np)
            spectrum["PEAKS_LIST"] = peak_list_np
            return spectrum
        else:
            spectrum['DELETION_REASON'] = "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number)."
            scripts.deletion_report.deleted_spectrum_list.append(spectrum)
            scripts.deletion_report.no_precursor_mz += 1
            return None
    elif "_GC" in spectrum["FILENAME"] or re.search(r"\bGC\b", spectrum["INSTRUMENTTYPE"]):
        float_precursor_mz = None
        peak_list_np = peak_list_cleaning(spectrum, peak_list, float_precursor_mz)

        if peak_list_np.size == 0:
            return None

        intensities = peak_list_np[:, 1]
        spectrum["ENTROPY"] = str(entropy_calculation(intensities))
        if parameters_dict["remove_spectrum_under_entropy_score"] == 1.0:
            if re.search(scripts.globals_vars.float_check_pattern, str(spectrum["ENTROPY"])):
                if float(spectrum["ENTROPY"]) < parameters_dict["remove_spectrum_under_entropy_score_value"]:
                    spectrum['DELETION_REASON'] = "spectrum deleted because it's entropy score is lower than the threshold selected by the user."
                    scripts.deletion_report.deleted_spectrum_list.append(spectrum)
                    scripts.deletion_report.low_entropy_score += 1
                    return None
        # If numpy array is empty, it returns none
        if peak_list_np.size == 0:
            return None
        spectrum["NUM PEAKS"] = str(peak_list_np.shape[0])
        # Convert numpy array back to string and update 'PEAKS_LIST' in spectrum
        peak_list_np = peak_list_to_str(peak_list_np)
        spectrum["PEAKS_LIST"] = peak_list_np
        return spectrum

    return spectrum


def write_deleted_spectrums(output_directory):
    """
    Groups the deleted spectrums by the DELETION_REASON and writes each group to a separate CSV file.

    This function creates a directory named 'DELETED_SPECTRUMS' inside the given output directory,
    if it does not already exist. It groups the spectrums by "DELETION_REASON" using if-elif for specific
    reasoning and then writes each group to a predefined filename set manually.

    Arguments:
    - output_directory (str): The path to the directory where the deleted spectrums directory
      and files should be created.

    Returns:
    - None
    """
    # Ensure the target directory exists
    deleted_spectrums_dir = os.path.join(output_directory, 'DELETED_SPECTRUMS')

    # Convert the list of dictionaries into a DataFrame
    spectrums_df = pd.DataFrame(scripts.deletion_report.deleted_spectrum_list)

    # Group the data by DELETION_REASON
    if "DELETION_REASON" in spectrums_df.columns:
        grouped = spectrums_df.groupby("DELETION_REASON")

        # Iterate through the groups explicitly checking for known reasons
        for reason, group in grouped:
            if reason == "spectrum deleted because peaks list is empty":
                file_name = "peaks_list_is_empty.csv"
            elif reason == "spectrum deleted because precursor mz is less than or equal to zero.":
                file_name = "precursor_mz_less_than_or_equal_zero.csv"
            elif reason == "spectrum deleted because it's entropy score is lower than the threshold selected by the user.":
                file_name = "entropy_score_lower_than_threshold.csv"
            elif reason == "spectrum deleted because precursor mz field is empty or contains invalid characters (not a floating number).":
                file_name = "precursor_mz_invalid_or_empty.csv"
            elif reason == "spectrum deleted because it has neither inchi nor smiles nor inchikey":
                file_name = "no_inchi_smiles_or_inchikey.csv"
            elif reason == "spectrum deleted because its adduct field is empty or the value entered is not an adduct":
                file_name = "adduct_empty_or_invalid.csv"
            elif reason == "spectrum deleted because the adduct corresponds to the wrong ionization mode (neg adduct in pos ionmode).":
                file_name = "wrong_adduct_neg_in_pos.csv"
            elif reason == "spectrum deleted because the adduct corresponds to the wrong ionization mode (pos adduct in neg ionmode).":
                file_name = "wrong_adduct_pos_in_neg.csv"
            elif reason == "spectrum deleted because its number of peaks is below the threshold chosen by the user":
                file_name = "number_of_peaks_below_threshold.csv"
            elif reason == "spectrum deleted because peaks list is empty after removing peaks above precursor m/z":
                file_name = "peaks_empty_after_above_precursor_mz_removal.csv"
            elif reason == "spectrum deleted because peaks list is empty after removing peaks out of mz range choiced by the use":
                file_name = "peaks_empty_after_mz_range_removal.csv"
            elif reason == "spectrum deleted because peaks list does not contain minimum number of high peaks required according to the value choiced by the user":
                file_name = "insufficient_high_peaks.csv"

            # Write the group to a CSV file
            file_path = os.path.join(deleted_spectrums_dir, file_name)
            group.to_csv(file_path, sep='\t', index=False, quotechar='"')

        # Clear the deleted spectrum list in memory
        scripts.deletion_report.deleted_spectrum_list = []


def spectrum_cleaning_processing(spectrum_list, output_directory, progress_callback=None, total_items_callback=None, prefix_callback=None,
                                 item_type_callback=None):
    """
    Main function used for performing spectrum cleaning operation on multiple spectrums.

    This function uses the concurrent.futures module to perform cleaning operation on multiple spectrums concurrently.
    The spectrum list is divided into chunks of determined size. Each chunk is then processed concurrently by utilizing
    the ThreadPoolExecutor.

    During the operation, progress is reported via callbacks instead of a progress bar. The callbacks are triggered
    as the cleaning process progresses.

    Finally, the cleaned spectrums are collected and returned.

    :param spectrum_list: A list of spectrums that need to be cleaned
    :type spectrum_list: list
    :param progress_callback: A function to update the progress (optional)
    :param total_items_callback: A function to set the total number of items (optional)
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional)
    :param item_type_callback: A function to specify the type of items being processed (optional)
    :return: A list of cleaned spectrums
    :rtype: list
    """

    # Define a task prefix via the callback (if provided)
    if prefix_callback:
        prefix_callback("cleaning spectrums:")

    # Specify the type of items being processed (if provided)
    if item_type_callback:
        item_type_callback("spectra")

    # Define the total number of spectrums for processing
    if total_items_callback:
        total_items_callback(len(spectrum_list), 0)  # Set total items and display initial completed = 0

    # Define the chunk size: the number of spectrums that can be processed at once
    chunk_size = calculate_maximized_chunk_size(data_list=spectrum_list)

    # Initialize the list that will hold the cleaned spectra
    final = []

    # Variable to track processed items
    processed_items = 0

    # Loop over the spectrum list, processing chunks of size `chunk_size`
    for i in range(0, len(spectrum_list), chunk_size):
        # Create a thread pool executor for concurrent processing
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Extract the current chunk
            chunk = spectrum_list[i:i + chunk_size]
            # Execute the `spectrum_cleaning` function for each item in the current chunk concurrently
            results = list(executor.map(spectrum_cleaning, chunk))

        # Collect all non-None results (cleaned spectrums)
        final.extend([res for res in results if res is not None])

        # Update the number of processed items
        processed_items += len(chunk)

        # Trigger the progress callback with the updated count
        if progress_callback:
            progress_callback(processed_items)

    write_deleted_spectrums(output_directory)

    # Return the final cleaned spectrum list
    return final
