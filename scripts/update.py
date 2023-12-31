import json
import re

def init_json_update_file(json_update_file):
    """
    Initializes a JSON update file if not provided.

    :param json_update_file: The existing JSON update file. If not provided, a new one will be created.
    :type json_update_file: dict or None
    :return: The initialized or existing JSON update file.
    :rtype: dict
    """
    if not json_update_file:
        json_update_file = {"FRAGHBID_LIST": []}
        return json_update_file
    else:
        return json_update_file


def check_for_update(spectrum_list):
    """
    :param spectrum_list: A list of spectrums to check for updates.
    :return: A list of spectrums that have a FRAGHUBID in the difference_list.
    """
    print("-- checking for updates --")
    # Ajout du code pour ouvrir et lire le fichier JSON
    with open('./datas/update.json', 'r') as f:
        json_update_file = json.load(f)

    json_update_file = init_json_update_file(json_update_file)

    # associer chaque spectrum avec son FRAGHUBID
    fraghub_id_spectrums = [(re.search("(?:FRAGHUBID: )(.*)", spectrum).group(1), spectrum) for spectrum in spectrum_list if re.search("(?:FRAGHUBID: )(.*)", spectrum)]
    fraghub_id_list = [fraghub_id for fraghub_id, spectrum in fraghub_id_spectrums]

    if not json_update_file:
        json_update_file = {"FRAGHBID_LIST": fraghub_id_list}
        # écrire les modifications dans le fichier JSON
        with open('./datas/update.json', 'w') as f:
            json.dump(json_update_file, f)
        return spectrum_list, False
    else:
        json_fraghub_id_list = json_update_file["FRAGHBID_LIST"]
        # récupérer la liste des id qui sont dans fraghub_id_list mais pas dans json_fraghub_id_list
        difference_list = list(set(fraghub_id_list) - set(json_fraghub_id_list))
        # ajouter les éléments de difference_list à json_fraghub_id_list
        json_fraghub_id_list.extend(difference_list)
        json_update_file["FRAGHBID_LIST"] = json_fraghub_id_list
        # écrire les modifications dans le fichier JSON
        with open('./datas/update.json', 'w') as f:
            json.dump(json_update_file, f)

        # retourner la liste des spectres qui ont un FRAGHUBID dans difference_list
        return [spectrum for fraghub_id, spectrum in fraghub_id_spectrums if fraghub_id in difference_list], True
