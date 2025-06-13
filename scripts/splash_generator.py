from scripts.calculate_maximized_chunk_size import *
from scripts.splash.spectrum_type import SpectrumType
from scripts.splash.spectrum import Spectrum
from scripts.splash.splash import Splash
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


def generate_splash_processing(spectrum_list, files, progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None):
    """
    Perform parallel processing of the given spectrum list and generate SPLASH for each spectrum,
    with support for progress reporting via callbacks.

    :param spectrum_list: A list of spectra to process.
    :param files: The name of the file related to spectrum_list (used for display purposes).
    :param progress_callback: A function to update progress (optional).
    :param total_items_callback: A function to set the total number of items (optional).
    :param prefix_callback: A function to dynamically set the prefix for the operation (optional).
    :param item_type_callback: A function to specify the type of items being processed (optional).
    :return: A list containing the SPLASH generated for each spectrum.
    """
    # Extract the file name from its path
    filename = os.path.basename(files).split("_")[0]

    # Set the prefix via the callback, if provided
    if prefix_callback:
        prefix_callback(f"generating SPLASH for [{filename}]:")

    # Specify the type of items via the callback
    if item_type_callback:
        item_type_callback("spectra")

    # Set the total via the callback (number of items in the spectrum list)
    if total_items_callback:
        total_items_callback(len(spectrum_list), 0)  # total = length of spectrum_list, completed = 0

    # Calculate the chunk size once
    chunk_size = calculate_maximized_chunk_size(data_list=spectrum_list)

    # List to store the final results
    final = []

    # Variable to track progress
    processed_items = 0

    # Divide the spectrum list into chunks and process each chunk
    for i in range(0, len(spectrum_list), chunk_size):
        # Create a chunk of spectra
        chunk = spectrum_list[i:i + chunk_size]

        # Use `ThreadPoolExecutor` for parallel processing
        with concurrent.futures.ThreadPoolExecutor() as executor:
            # Apply the `generate_splash` function to each spectrum in the chunk
            results = list(executor.map(generate_splash, chunk))

        # Filter results to exclude `None` values
        final.extend([res for res in results if res is not None])

        # Update the number of processed items
        processed_items += len(chunk)

        # Update progress via the callback, if provided
        if progress_callback:
            progress_callback(processed_items)

    # Return the final list of results
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
