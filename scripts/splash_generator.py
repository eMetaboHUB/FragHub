from splash import Spectrum, SpectrumType, Splash
from calculate_maximized_chunk_size import *
import concurrent.futures
import os

def hash_spectrum_data(spectrum_data):
    """
    Hashes the spectrum data with SPLASH key.

    :param spectrum_data: The spectrum data.
    :return: The hashed spectrum.
    """
    peak_list = [tuple(peaks) for peaks in spectrum_data["PEAKS_LIST"]]

    # Check if peak_list exist
    if peak_list:
        # calculate splash key
        try:
            spectrum = Spectrum(peak_list, SpectrumType.MS)
            return Splash().splash(spectrum)
        except:
            return None
    else:
        return None

def generate_splash(spectrum):
    """

    :param spectrum: a dictionary containing the spectrum data
    :return: a dictionary containing the spectrum data updated with a new splash

    This method generates a splash for the given spectrum data. The spectrum data is hashed using the `hash_spectrum_data` method and the resultant hash is converted into a string.

    If no splash can be generated (i.e., the hash is empty), None is returned.

    The generated splash is then added to the spectrum data dictionary using the key "SPLASH".

    Finally, the spectrum data with the newly added splash is returned.

    """
    # Hash the spectrum data and convert the resultant hash into a string.
    # splash holds the hashed id of the spectrum data
    if not isinstance(spectrum, dict):
        return None

    splash = str(hash_spectrum_data(spectrum))

    # Return None if splash is empty i.e., no splash could be generated.
    if not splash:
        return None

    # Add the generated splash to the spectrum data dictionary.
    # This line is updating the spectrum dictionary with the newly generated splash
    spectrum["SPLASH"] = splash

    # Return the spectrum data with the newly added splash
    # This will be a dictionary containing the spectrum data, updated with a new splash.
    return spectrum


def generate_splash_processing(spectrum_list, files, progress_callback=None, total_items_callback=None,
                               prefix_callback=None, item_type_callback=None):
    """
    Perform parallel processing of the given spectrum list and generate splash for each spectrum,
    with support for progress reporting via callbacks.

    :param spectrum_list: A list of spectra to process.
    :param files: The name of the file related to spectrum_list (used for display purposes).
    :param progress_callback: A function to update the progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items processed (optional).
    :return: A list containing the SPLASH generated for each spectrum.
    """
    # Extraire le nom du fichier depuis son chemin
    filename = os.path.basename(files).split("_")[0]

    # Définir le préfixe via le callback, si fourni
    if prefix_callback:
        prefix_callback(f"generating SPLASH for [{filename}]:")

    # Spécifier le type d'éléments via le callback
    if item_type_callback:
        item_type_callback("spectra")

    # Définir le total via le callback (nombre d'éléments dans la liste des spectres)
    if total_items_callback:
        total_items_callback(len(spectrum_list), 0)  # total = longueur de spectrum_list, completed = 0

    # Taille des chunks pour le traitement par lots
    chunk_size = calculate_maximized_chunk_size(data_list=spectrum_list)

    # Liste pour stocker les résultats finaux
    final = []

    # Variable pour suivre la progression
    processed_items = 0

    # Diviser la liste des spectres en chunks et traiter chaque chunk
    for i in range(0, len(spectrum_list), chunk_size):
        # Créer un chunk de spectres
        chunk = spectrum_list[i:i + chunk_size]

        # Utiliser `ThreadPoolExecutor` pour le traitement parallèle
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Appliquer la fonction `generate_splash` à chaque spectre dans le chunk
            results = list(executor.map(generate_splash, chunk))

        # Filtrer les résultats pour exclure les valeurs `None`
        final.extend([res for res in results if res is not None])

        # Mettre à jour le nombre d'éléments traités
        processed_items += len(chunk)

        # Mettre à jour la progression via le callback, si fourni
        if progress_callback:
            progress_callback(processed_items)

    # Retourner la liste finale des résultats
    return final


# process_converted_after function processes the converted spectrum list.
# Based on the mode, it identifies the appropriate file path for the converted file.
def process_converted_after(spectrum_list, mode, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Process converted spectrum list after conversion.
    :param spectrum_list: List of converted spectra.
    :param mode: Conversion mode. Can be "MSP", "XML", "CSV", or "JSON".
    :return: Processed spectrum list.
    """
    # Assign relevant file path based on conversion mode
    if mode == "MSP":
        file_path = os.path.abspath("../INPUT/CONVERTED/MSP_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "XML":
        file_path = os.path.abspath("../INPUT/CONVERTED/XML_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "CSV":
        file_path = os.path.abspath("../INPUT/CONVERTED/CSV_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "JSON":
        file_path = os.path.abspath("../INPUT/CONVERTED/JSON_converted.json")
        filename = os.path.basename(file_path)
    elif mode == "MGF":
        file_path = os.path.abspath("../INPUT/CONVERTED/MGF_converted.json")
        filename = os.path.basename(file_path)

    # Generate splash for the spectrum list
    spectrum_list = generate_splash_processing(spectrum_list, filename, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Return the processed spectrum list
    return spectrum_list

def generate_splash_id(FINAL_MSP, FINAL_CSV, FINAL_JSON, FINAL_MGF, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Process the converted files and generate a splash ID for each file type.
    - FINAL_MSP: The path to the converted MSP file.
    - FINAL_XML: The path to the converted XML file.
    - FINAL_CSV: The path to the converted CSV file.
    - FINAL_JSON: The path to the converted JSON file.
    - FINAL_MGF: The path to the converted MGF file.

    - Return: A tuple containing the splash IDs for each file type (MSP, XML, CSV, JSON, MGF).
    """

    # Check if the MSP file path is valid, if it is, process the file to generate splash id
    if FINAL_MSP:
        FINAL_MSP = process_converted_after(FINAL_MSP, "MSP", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Similar check for CSV file path
    if FINAL_CSV:
        FINAL_CSV = process_converted_after(FINAL_CSV, "CSV", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Similar check for JSON file path
    if FINAL_JSON:
        FINAL_JSON = process_converted_after(FINAL_JSON, "JSON", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Similar check for MGF file path
    if FINAL_MGF:
        FINAL_MGF = process_converted_after(FINAL_MGF, "MGF", progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

    # Finally return a tuple containing splash ids for each file type
    return FINAL_MSP, FINAL_CSV, FINAL_JSON, FINAL_MGF
