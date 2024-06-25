from splash import Spectrum, SpectrumType, Splash
import concurrent.futures
from tqdm import tqdm
import os

def hash_spectrum_data(spectrum_data):
    """
    Hashes the spectrum data with SPLASH key.

    :param spectrum_data: The spectrum data.
    :return: The hashed spectrum.
    """
    # Convert spectrum data to string

    # Search for inchikey_update_pattern in spectrum data
    peak_list = [tuple(peaks) for peaks in spectrum_data["PEAKS_LIST"]]

    # Check if both inchikey and peak_list exist
    if peak_list:
        # Combine inchikey and peak list into one string with a newline separator
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

def generate_splash_processing(spectrum_list, files):
    """
    Perform parallel processing of the given spectrum list and generate splash for each spectrum.
    :param spectrum_list: A list of spectra.
    :return: A list of splash generated for each spectrum.
    """
    filename = os.path.basename(files)  # Extract the base name of the file path
    chunk_size = 5000  # Set the size of chunks
    final = []  # Create an empty list to store final results
    # Create a progress bar with the total length of spectrum_list and a relevant description
    progress_bar = tqdm(total=len(spectrum_list), unit=" spectrums", colour="green", desc="{:>70}".format(f"generating SPLASH on [{filename}]"))
    # Dividing the spectrum list into chunks
    for i in range(0, len(spectrum_list), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = spectrum_list[i:i + chunk_size]  # create a chunk
            # Use a executor's map function to apply 'generate_splash' function to each spectrum in the chunk, and convert it to a list
            results = list(executor.map(generate_splash, chunk))
            progress_bar.update(len(chunk))  # update the progress bar by the size of the chunk processed
        final.extend([res for res in results if res is not None])  # Extend the final results list with results that are not None
    progress_bar.close()  # close the progress bar
    return final  # return the final results

# process_converted_after function processes the converted spectrum list.
# Based on the mode, it identifies the appropriate file path for the converted file.
def process_converted_after(spectrum_list, mode):
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
    spectrum_list = generate_splash_processing(spectrum_list, filename)

    # Return the processed spectrum list
    return spectrum_list

def generate_splash_id(FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF):
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
        FINAL_MSP = process_converted_after(FINAL_MSP, "MSP")

    # Similar check for XML file path
    if FINAL_XML:
        FINAL_XML = process_converted_after(FINAL_XML, "XML")

    # Similar check for CSV file path
    if FINAL_CSV:
        FINAL_CSV = process_converted_after(FINAL_CSV, "CSV")

    # Similar check for JSON file path
    if FINAL_JSON:
        FINAL_JSON = process_converted_after(FINAL_JSON, "JSON")

    # Similar check for MGF file path
    if FINAL_MGF:
        FINAL_MGF = process_converted_after(FINAL_MGF, "MGF")

    # Finally return a tuple containing splash ids for each file type
    return FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF
