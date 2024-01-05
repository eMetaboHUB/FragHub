import concurrent.futures
from tqdm import tqdm
import json
import re

global fraghub_id_pattern
fraghub_id_pattern = re.compile("(?:FRAGHUBID: )(.*)")



def init_json_update_file(json_update_file):
    """
    :param json_update_file: The JSON update file to be initialized. It should be a dictionary or None.
    :return: The initialized JSON update file. If the input JSON update file is None, a new dictionary will be created with the key "FRAGHBID_LIST" and an empty list as its value. If the
    * input JSON update file is not None, it will be returned as is.
    """
    first_run = False

    if not bool(json_update_file):
        first_run = True
        json_update_file = {"FRAGHBID_LIST": {}}
        return json_update_file, first_run
    else:
        return json_update_file, first_run

def check_for_update(spectrum):
    """
    :param spectrum: The spectrum to check for updates.
    :return: If the spectrum is not found in the dictionary or the dictionary is empty, returns the spectrum and the fraghub_id_spectrum. Otherwise, returns None.
    """
    fraghub_id_spectrum = re.search(fraghub_id_pattern, spectrum)
    if fraghub_id_spectrum:
        fraghub_id_spectrum = fraghub_id_spectrum.group(1)

    fraghub_id_dict = json_update_file["FRAGHBID_LIST"]

    # Vérifier si le dictionnaire est vide ou si fraghub_id_spectrum n'est pas dans le dictionnaire.
    if not fraghub_id_dict or fraghub_id_spectrum not in fraghub_id_dict:
        return spectrum, fraghub_id_spectrum
    else:
        return None

def check_for_update_processing(spectrum_list):
    """

    :param spectrum_list: A list of spectrums to check for updates.
    :return: A list of spectrums that have been updated.

    """
    global json_update_file

    with open('../datas/update.json', 'r') as f:
        json_update_file = json.load(f)

    json_update_file, first_run = init_json_update_file(json_update_file)

    chunk_size = 5000
    final = []
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green",
                        desc="{:>80}".format("checking for updates"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(spectrum_list), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = spectrum_list[i:i + chunk_size]
            results = list(executor.map(check_for_update, chunk))
            progress_bar.update(len(chunk))

        final.extend([res for res in results if res is not None])

    final_spectrum_list = [res[0] for res in final]
    new_fraghubid = {res[1]: True for res in final}

    if final:
        update = True
    else:
        update = False

    json_update_file["FRAGHBID_LIST"].update(new_fraghubid)

    # écrire les modifications dans le fichier JSON
    with open('../datas/update.json', 'w') as f:
        json.dump(json_update_file, f, ensure_ascii=False, indent=4)

    progress_bar.close()

    return final_spectrum_list, update, first_run
