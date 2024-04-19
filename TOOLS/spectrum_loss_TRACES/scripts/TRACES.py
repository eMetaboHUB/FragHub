import concurrent.futures
from tqdm import tqdm
import pandas as pd

def structure_traces_one(spectrums):
    """
    :param spectrums: A list containing values related to spectra
    :return: A dictionary containing processed values from the input list
    """

    if spectrums[0] is None:
        temp_spectrums = spectrums[5]
        temp_spectrums.update({"no_smiles_no_inchi": spectrums[1],
                               "no_or_bad_precursormz_and_no_or_bad_addcut": spectrums[2],
                               "no_peaks_list": spectrums[3],
                               "peaks_filters_not_valid": spectrums[4],
                               "will_be_deleted_because_of_RDkit": False
                               })

        # Used list-comprehension to ensure every value is a list
        # temp_spectrums = {k: [v] for k, v in temp_spectrums.items()}

        return temp_spectrums

def spectrum_traces(spectrum_list, profile_name):
    """
    Main function used for performing spectrum cleaning operation on multiple spectrums.

    This function uses the concurrent.futures module to perform cleaning operation on multiple spectrums concurrently.
    The spectrum list is divided into chunks of determined size. Each chunk is then processed concurrently by utilizing
    the ThreadPoolExecutor.

    During the operation, a progress bar is shown to indicate the progress of the operation. The progress bar is updated
    every time a chunk of spectrums has been processed.

    Finally, the cleaned spectrums are collected and returned.

    :param spectrum_list: A list of spectrums that need to be cleaned
    :type spectrum_list: list
    :return: A list of cleaned spectrums
    :rtype: list
    """

    # Define the chunk size: the number of spectrums that can be processed at once
    chunk_size = 5000
    # Initialize the list that will hold the cleaned spectra
    final = []
    # Initialize a progress bar to monitor the processing
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>70}".format("structure trackers"))
    # Loop over the spectrum list using the defined chunk size
    for i in range(0, len(spectrum_list), chunk_size):
        # Create a ThreadPoolExecutor for concurrent processing
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Extract a chunk of spectrum from the spectrum list
            chunk = spectrum_list[i:i + chunk_size]
            # Execute the spectrum_cleaning function on each spectrum from the chunk concurrently
            results = list(executor.map(structure_traces_one, chunk))
            # Update the progress bar
            progress_bar.update(len(chunk))
        # Collect all non-None cleaned spectrums from the chunk and append them to the final cleaned spectrum list
        final.extend([res for res in results if res is not None])
    # Close the progress bar after all has been processed
    progress_bar.close()
    return final  # Return the cleaned spectrum list

def structure_traces_two(spectrum_list_TRACES_DF, DELETED_CONCATENATE_DF, profile_name):
    """
    :param spectrum_list_TRACES: A list containing dictionaries of spectrums. Each dictionary represents a spectrum and should have keys 'old_spectrum', 'no_smiles_no_inchi', 'no_or_bad_precursormz_and_no_or_bad_addcut', 'no_peaks_list', and 'minimum_peaks_not_requiered'. The value of 'old_spectrum' key should be a dictionary containing spectrum information.
    :param DELETED_CONCATENATE_DF: A DataFrame containing deleted spectrums.
    :param profile_name: A string representing the profile name.
    :return: None

    This method takes in a list of spectrums, merges them into a DataFrame, concatenates it with DELETED_CONCATENATE_DF, and saves the resulting DataFrame as an Excel file in the specified directory.
    """
    DF = pd.concat([spectrum_list_TRACES_DF, DELETED_CONCATENATE_DF])

    DF.to_excel(rf"../OUTPUT/{profile_name}/TRACKER_2.xlsx", index=False)


