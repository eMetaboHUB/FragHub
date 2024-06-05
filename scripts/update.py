import concurrent.futures
from tqdm import tqdm
import json

def init_json_update_file(json_update_file):
    """
    :param json_update_file: The JSON update file to be initialized. It should be a dictionary or None.
    :return: The initialized JSON update file. If the input JSON update file is None, a new dictionary will be created with the key "SPLASH_LIST" and an empty list as its value. If the
    * input JSON update file is not None, it will be returned as is.
    """
    # Initialize the "first run" flag as False
    first_run = False

    # If the json_update_file is not provided (None), set the "first run" flag to True,
    # and initialize a new dictionary with "SPLASH_LIST" key and an empty value,
    # then return this dictionary and the first_run flag
    if not bool(json_update_file):
        first_run = True
        # Creating a new dictionary to represent the JSON update file
        json_update_file = {"SPLASH_LIST": {}}
        return json_update_file, first_run
    else:
        # If json_update_file is already provided, simply return it as is along with the first_run flag
        return json_update_file, first_run

def check_for_update(spectrum):
    """
    :param spectrum: The spectrum to check for updates.
    :return: If the spectrum is not found in the dictionary or the dictionary is empty, returns the spectrum and the fraghub_id_spectrum. Otherwise, returns None.
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
        # If the SPLASH is found in the dictionary, an update is not needed. Return None
        return None

def check_for_update_processing(spectrum_list, profile_name):
    """
    This method checks for updates in the given spectrum list using a profile name.
    :param spectrum_list: A list of spectrums to check for updates.
    :param profile_name: The name of the profile to use for checking updates.
    :return: A tuple containing the updated spectrum list, a flag indicating if an update was found, and the first run flag.
    """
    global json_update_file
    # Read json file that contains previous update status
    with open(f'../datas/updates/{profile_name}.json', 'r') as f:
        json_update_file = json.load(f)

    # Initialize json file and get first run flag
    json_update_file, first_run = init_json_update_file(json_update_file)

    # chunk size for how many spectrums will be checked simultaneously
    chunk_size = 5000

    final = []  # final list of processed spectrums

    # Progress bar for visual representation of the task
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>70}".format("checking for updates"))

    # Iterating the spectrums in chunk size
    for i in range(0, len(spectrum_list), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = spectrum_list[i:i + chunk_size]

            # Checking each spectrum chunk in concurrent manner
            results = list(executor.map(check_for_update, chunk))
            progress_bar.update(len(chunk))

        # Adding successfully checked spectrums to the final list
        final.extend([res for res in results if res is not None])

    # Creating new list with updated spectrums only
    final_spectrum_list = [res[0] for res in final]

    # Newly found SPLASH list
    new_splash = {res[1]: True for res in final}

    # Update flag - if new spectrums found, set it as True, otherwise, False
    if final:
        update = True
    else:
        update = False

    # Updating json file with new SPLASH list
    json_update_file["SPLASH_LIST"].update(new_splash)

    # Write updated data to the json file
    with open(f'../datas/updates/{profile_name}.json', 'w') as f:
        json.dump(json_update_file, f, ensure_ascii=False, indent=4)

    # Close the progress bar after finishing the task
    progress_bar.close()

    # Return results : updated spectrum list, update flag, first run flag
    return final_spectrum_list, update, first_run
