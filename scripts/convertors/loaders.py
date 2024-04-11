from tqdm import tqdm
import ijson
import os
import re

def load_spectrum_list_from_msp(msp_file_path):
    """
    Load a spectrum list from a given MSP (Mass Spectral Peak) file.

    :param msp_file_path: The path to the MSP file.
    :return: The list of spectra read from the file. Each spectrum is represented as a string.

    """
    filename = os.path.basename(msp_file_path)
    spectrum_list = []
    buffer = [f"FILENAME: {os.path.basename(msp_file_path)}\n"]  # initialisation of buffer with the filename

    total_size = os.path.getsize(msp_file_path)  # get the total size of the file in bytes

    with open(msp_file_path, 'r', encoding="UTF-8") as file:
        for line in tqdm(file, total=total_size, unit="B", unit_scale=True, colour="green", desc="{:>70}".format(f"loading [{filename}]")): # wrap this with tqdm
            if line.strip() == '':
                if buffer:
                    spectrum = '\n'.join(buffer)
                    spectrum = re.sub(r"FILENAME: .*\n", f"FILENAME: {os.path.basename(msp_file_path)}\n", spectrum, flags=re.IGNORECASE)
                    spectrum_list.append(spectrum)
                    buffer = [f"FILENAME: {os.path.basename(msp_file_path)}\n"]  # buffer reinitialisation with the filename
            else:
                buffer.append(line.strip())

    return spectrum_list

def load_spectrum_list_from_mgf(mgf_file_path):
    """
    Load a spectrum list from a given MGF (Mascot Generic Format) file.

    :param mgf_file_path: The path to the MGF file.
    :return: The list of spectra read from the file. Each spectrum is represented as a string.

    """
    filename = os.path.basename(mgf_file_path)
    spectrum_list = []
    buffer = [f"FILENAME={os.path.basename(mgf_file_path)}"]  # initialise the buffer with the filename

    total_size = os.path.getsize(mgf_file_path)  # get the total size of the file in bytes

    with open(mgf_file_path, 'r', encoding="UTF-8") as file:
        for line in tqdm(file, total=total_size, unit="B", unit_scale=True, colour="green", desc="{:>70}".format(f"loading [{filename}]")):  # wrap this with tqdm
            if line.strip() == 'END IONS':
                if buffer:
                    spectrum = '\n'.join(buffer)
                    spectrum = re.sub(r"FILENAME=.*\n", f"FILENAME={os.path.basename(mgf_file_path)}\n", spectrum, flags=re.IGNORECASE)  # remove any existing 'FILENAME=' line
                    spectrum_list.append(spectrum)
                    buffer = [f"FILENAME={os.path.basename(mgf_file_path)}"]  # reinitialise the buffer with the filename
            else:
                buffer.append(line.strip())

    return spectrum_list

def load_spectrum_list_json(json_file_path):
    """
    Load spectra from a JSON file and return a generator of spectra.

    :param json_file_path: Path to the JSON file.
    :return: Generator of spectra.
    """

    # First, calculate total bytes for tqdm
    filename = os.path.basename(json_file_path)
    total_bytes = os.path.getsize(json_file_path)

    with open(json_file_path, 'r', encoding="UTF-8") as file:
        # ijson.items(file, 'item') returns a generator yielding items in a JSON file
        spectra = ijson.items(file, 'item')
        # Create tqdm progress bar
        progress = tqdm(total=total_bytes, unit="B", unit_scale=True, colour="green", desc="{:>70}".format(f"loading [{filename}]"))

        filename = os.path.basename(json_file_path)
        for spectrum in spectra:
            spectrum["filename"] = filename
            progress.update(len(str(spectrum)))
            yield spectrum

        progress.close()