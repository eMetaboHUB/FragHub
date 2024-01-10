from tqdm import tqdm
import ijson
import os

def load_spectrum_list_from_msp(msp_file_path):
    """
    Load a spectrum list from a given MSP (Mass Spectral Peak) file.

    :param msp_file_path: The path to the MSP file.
    :return: The list of spectra read from the file. Each spectrum is represented as a string.

    Example usage:
    ```
    msp_file_path = "path/to/spectrum.msp"
    spectrum_list = load_spectrum_list(msp_file_path)
    print(spectrum_list)
    ```
    """
    filename = os.path.basename(msp_file_path)
    spectrum_list = []
    buffer = []

    total_size = os.path.getsize(msp_file_path)  # get the total size of the file in bytes

    with open(msp_file_path, 'r', encoding="UTF-8") as file:
        for line in tqdm(file, total=total_size, unit="B", unit_scale=True, colour="green", desc="{:>80}".format(f"loading [{filename}]")): # wrap this with tqdm
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

def load_spectrum_list_json(json_file_path):
    """
    Load spectra from a JSON file and return a list of spectra.

    :param json_file_path: Path to the JSON file.
    :return: List of spectra.
    """
    spectrum_list = []

    # First, calculate total bytes for tqdm
    total_bytes = os.path.getsize(json_file_path)

    with open(json_file_path, 'r', encoding="UTF-8") as file:
        # ijson.items(file, 'item') returns a generator yielding items in a JSON file
        spectra = ijson.items(file, 'item')
        # Create tqdm progress bar
        progress = tqdm(total=total_bytes, unit="B", unit_scale=True, colour="green", desc="{:>80}".format("Loading file"))

        for spectrum in spectra:
            spectrum_list.append(spectrum)
            progress.update(len(str(spectrum)))
        progress.close()

    return spectrum_list