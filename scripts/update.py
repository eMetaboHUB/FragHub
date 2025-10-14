from scripts.calculate_maximized_chunk_size import *
import concurrent.futures
import scripts.deletion_report
import pandas as pd
import os.path
import json

def init_json_update_file(json_update_file):
    """
    :param json_update_file: The JSON update file to be initialized. It should be a dictionary or None.
    :return: The initialized JSON update file. If the input JSON update file is None, a new dictionary will be created with the key "SPLASH_LIST" and an empty list as its value. If the
    * input JSON update file is not None, it will be returned as is.
    """
    # Initialize the "first run" flag as False

    # If the json_update_file is not provided (None), set the "first run" flag to True,
    # and initialize a new dictionary with "SPLASH_LIST" key and an empty value,
    # then return this dictionary and the first_run flag
    if not bool(json_update_file):
        # Creating a new dictionary to represent the JSON update file
        json_update_file = {"SPLASH_LIST": {}}
        return json_update_file
    else:
        # If json_update_file is already provided, simply return it as is along with the first_run flag
        return json_update_file



def check_for_update(spectrum):
    """
    :param spectrum: The spectrum to check for updates.
    :return: If the spectrum is not found in the dictionary or the dictionary is empty, returns the spectrum and the fraghub_id_spectrum. Otherwise, appends the spectrum to the deleted_spectrum_list and returns None.
    """
    # Extract the SPLASH from the provided spectrum
    fraghub_id_spectrum = spectrum["SPLASH"]

    # Get the dictionary of SPLASH from the json file
    fraghub_id_dict = json_update_file["SPLASH_LIST"]

    # Check if the fraghub_id_dict is empty or if fraghub_id_spectrum does not exist in the fraghub_id_dict
    if not fraghub_id_dict or fraghub_id_spectrum not in fraghub_id_dict:
        # If the dictionary is empty, or the SPLASH is not present in the dictionary, an update is needed. Return the spectrum and the SPLASH
        return spectrum, fraghub_id_spectrum
    else:
        spectrum['DELETION_REASON'] = "spectrum deleted because already processed in a previous run."
        # If the SPLASH is found in the dictionary, an update is not needed.
        # Append the spectrum to the deleted_spectrum_list
        scripts.deletion_report.deleted_spectrum_list.append(spectrum)
        # Return None since no update is needed
        return None


def check_for_update_processing(spectrum_list, output_directory, progress_callback=None, total_items_callback=None,
                                prefix_callback=None, item_type_callback=None):
    """
    Check for updates in the given spectrum list using a profile name, with progress reporting via callbacks.

    :param spectrum_list: A list of spectrums to check for updates.
    :param profile_name: The name of the profile to use for checking updates.
    :param progress_callback: A function to update the progress (processed items).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items (optional).
    :return: A tuple (updated spectrum list, update flag, first run flag).
    """
    global json_update_file

    total = len(spectrum_list)

    # Extra task prefix
    if prefix_callback:
        prefix_callback("checking for updates:")

    # Specify the type of item being processed
    if item_type_callback:
        item_type_callback("spectra")

    # Define the total number of items using the callback
    if total_items_callback:
        total_items_callback(len(spectrum_list), 0)  # Total = len(spectrum_list), Completed = 0

    update_file_path = os.path.join(output_directory, "updates.json")
    # Open the JSON file containing the previous update status
    with open(update_file_path, 'r') as f:
        json_update_file = json.load(f)

    # Initialize the update file and first run flag
    json_update_file = init_json_update_file(json_update_file)

    # Define chunk size for parallel processing
    chunk_size = calculate_maximized_chunk_size(data_list=spectrum_list)

    # Final list to store successfully checked spectra
    final = []
    processed_items = 0  # Track the number of processed items

    # Iterate through the spectrum list in chunks
    for i in range(0, len(spectrum_list), chunk_size):
        chunk = spectrum_list[i:i + chunk_size]

        # Use ThreadPoolExecutor for parallel processing
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Check each spectrum in the chunk
            results = list(executor.map(check_for_update, chunk))

        # Add successfully checked spectra to the final list
        final.extend([res for res in results if res is not None])

        # Update the number of processed items
        processed_items += len(chunk)

        # Update progress callback if provided
        if progress_callback:
            progress_callback(processed_items)

    # Create a cleaned list with only updated spectra
    final_spectrum_list = [res[0] for res in final]

    # Generate a new SPLASH list
    new_splash = {res[1]: True for res in final}

    # Determine whether an update occurred
    update = bool(final)

    # Update the SPLASH list in the JSON file
    json_update_file["SPLASH_LIST"].update(new_splash)

    # Write the updated JSON file back to disk
    with open(update_file_path, 'w') as f:
        json.dump(json_update_file, f, ensure_ascii=False, indent=4)

    scripts.deletion_report.previously_cleaned = total - len(final_spectrum_list)

    deleted_spectrums_dir = os.path.join(output_directory, 'DELETED_SPECTRUMS')
    previously_cleaned_file = os.path.join(deleted_spectrums_dir, 'previously_cleaned.csv')
    deleted_spectra_df = pd.DataFrame(scripts.deletion_report.deleted_spectrum_list)
    deleted_spectra_df.to_csv(previously_cleaned_file, sep='\t', index=False, quotechar='"')

    # Reinitialize the deleted_spectrum_list to free memory
    scripts.deletion_report.deleted_spectrum_list = []

    # Return the final spectrum list, the update flag, and the first run flag
    return final_spectrum_list, update
