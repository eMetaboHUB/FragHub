import concurrent.futures
from tqdm import tqdm
import hashlib
import os

def load_spectrum_list(msp_file_path):
    """
    Load spectra from an MSP file and return a list of spectra.

    :param msp_file_path: Path to the MSP file.
    :return: List of spectra.
    """
    spectrum_list = []
    buffer = []

    total_lines = sum(1 for line in open(msp_file_path, 'r', encoding="UTF-8")) # count the total number of lines in the file

    with open(msp_file_path, 'r', encoding="UTF-8") as file:
        for line in tqdm(file, total=total_lines, unit=" rows", colour="green", desc="\t     reading"): # wrap this with tqdm
            if line.strip() == '':
                if buffer:
                    spectrum_list.append('\n'.join(buffer))
                    buffer = []
            else:
                if not buffer:
                    buffer.append(f"FILENAME: {os.path.basename(msp_file_path)}") # adding filename to spectrum
                buffer.append(line.strip())

    # Add the last spectrum to the list
    if buffer:
        spectrum_list.append('\n'.join(buffer))

    return spectrum_list

def hash_spectrum_data(spectrum_data):
    """
    :param spectrum_data: The spectrum data to be hashed.
    :return: The SHA-256 hash of the spectrum data as a hexadecimal string.

    This method takes the spectrum data and converts it into a string. It then creates a SHA-256 object and updates it with the spectrum string encoded in UTF-8. Finally, it returns the
    * hexadecimal representation of the SHA-256 hash.
    """
    # Convertir le spectre data en une chaîne de caractères
    spectrum_string = str(spectrum_data)

    # Créer un objet sha256
    sha256 = hashlib.sha256()

    # Fournir les données de spectre à sha256
    sha256.update(spectrum_string.encode('utf-8'))

    # Retourner le hash sha256 en hex
    return sha256.hexdigest()

def genrate_fraghubid(spectrum):
    """
    Generate FragHubID for a given spectrum.

    :param spectrum: The spectrum data.
    :return: The spectrum data with Fragment Hub ID.

    Example:
    --------
    >>> spectrum = "MS/MS Spectrum"
    >>> genrate_fraghubid(spectrum)
    'FRAGHUBID: 123456\nMS/MS Spectrum'
    """
    hash_key = hash_spectrum_data(spectrum)
    spectrum = f"FRAGHUBID: {str(hash_key)}\n" + spectrum

    return spectrum
def generate_fraghub_id(msp_directory_path):
    """
    Generate the FragHub ID for each spectrum in the given MSP directory path.

    :param msp_directory_path: The path to the MSP directory.
    :type msp_directory_path: str
    :return: None
    """
    for files in os.listdir(msp_directory_path):
        if files.endswith(".msp"):
            msp_file_path = os.path.join(msp_directory_path, files)
            spectrum_list= load_spectrum_list(msp_file_path)
            spectrum_list = genrate_fraghubid_processing(spectrum_list)

def genrate_fraghubid_processing(spectrum_list):
    """
    Perform parallel processing of the given spectrum list and generate fraghubid for each spectrum.

    :param spectrum_list: A list of spectra.
    :return: A list of fraghubids generated for each spectrum.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(genrate_fraghubid, spectrum_list), total=len(spectrum_list), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.
