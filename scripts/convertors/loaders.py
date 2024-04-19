from tqdm import tqdm
import ijson
import os
import re

def load_spectrum_list_from_msp(msp_file_path):
    """
    Load a spectrum list from a given MSP (Mass Spectral Peak) file.

    Arguments:
    msp_file_path (str): The path to the MSP file.

    Returns:
    spectrum_list (List[str]): The list of spectra read from the file. Each spectrum is represented as a string.

    """

    # Count number of spectra
    num_spectra = sum(1 for line in open(msp_file_path, 'r', encoding="UTF-8") if line.strip() == '')

    # Get the name of the file
    filename = os.path.basename(msp_file_path)

    # Initialize a list to store spectra
    spectrum_list = []

    # Initialize a buffer to temporarily hold data
    buffer = [f"FILENAME: {os.path.basename(msp_file_path)}\n"]  # initialization of buffer with the filename

    # Open the file for reading
    with tqdm(total=num_spectra, unit="spectra", colour="green", desc="{:>70}".format(f"loading [{filename}]")) as pbar:
        with open(msp_file_path, 'r', encoding="UTF-8") as file:

            # For each line in the file
            for line in file:

                # If the line is empty
                if line.strip() == '':

                    # Check if buffer is not empty
                    if buffer:
                        # Join the data in the buffer into a single string, and add it to the list of spectra
                        spectrum = '\n'.join(buffer)
                        spectrum = re.sub(r"FILENAME: .*\n", f"FILENAME: {os.path.basename(msp_file_path)}\n", spectrum, flags=re.IGNORECASE)
                        spectrum_list.append(spectrum)

                        # Reset the buffer
                        buffer = [f"FILENAME: {os.path.basename(msp_file_path)}\n"]

                        # update the progress bar
                        pbar.update()

                # If the line is not empty, add it to the buffer.
                else:
                    buffer.append(line.strip())

    # Return the list of spectra
    return spectrum_list

def load_spectrum_list_from_mgf(mgf_file_path):
    """
    Load a spectrum list from a given MGF (Mascot Generic Format) file.
    :param mgf_file_path: The path to the MGF file.
    :return: The list of spectra read from the file. Each spectrum is represented as a string.
    """
    # Count the total number of spectra
    num_spectra = sum(1 for line in open(mgf_file_path, 'r', encoding="UTF-8") if line.strip() == 'END IONS')

    filename = os.path.basename(mgf_file_path)
    spectrum_list = []
    buffer = [f"FILENAME={filename}"]

    with tqdm(total=num_spectra, unit="spectra", colour="green", desc="{:>70}".format(f"loading [{filename}]")) as pbar:
        with open(mgf_file_path, 'r', encoding="UTF-8") as file:
            for line in file:
                if line.strip() == 'END IONS':
                    if buffer:
                        spectrum = '\n'.join(buffer)
                        spectrum = re.sub(r"FILENAME=.*\n", f"FILENAME={filename}\n", spectrum, flags=re.IGNORECASE)
                        spectrum_list.append(spectrum)
                        buffer = [f"FILENAME={filename}"]
                    pbar.update(1)  # Update progress bar each time a spectrum is being processed
                else:
                    buffer.append(line.strip())

    # return the list of spectra
    return spectrum_list

def load_spectrum_list_json(json_file_path):
    """
    This function is used to load a list of spectra from a given JSON file.

    :param json_file_path: A string representing the path to the JSON file containing the spectra.

    :return: A generator yielding each spectrum (a dictionary) from the JSON file, one at a time.
    """

    # Extract the filename from the given path
    filename = os.path.basename(json_file_path)

    # Calculate the total size of the file in bytes for progress reporting
    total_bytes = os.path.getsize(json_file_path)

    # Open the JSON file for reading
    # Note: The encoding is set to UTF-8 to support wide range of characters.
    with open(json_file_path, 'r', encoding="UTF-8") as file:
        # Create a generator to yield each 'item' in the JSON file using ijson.
        # ijson will yield each 'item' as a separate dictionary.
        spectra = ijson.items(file, 'item')

        # Initialize a progress bar for tracking file loading progress with tqdm.
        # The total size is set to total_bytes (file size), and a description is attached showing the filename being loaded.
        progress = tqdm(total=total_bytes, unit="B", unit_scale=True, colour="green", desc="{:>70}".format(f"loading [{filename}]"))

        # Loop over each spectrum in the file
        for spectrum in spectra:
            # Add the originating filename to the spectrum dictionary
            spectrum["filename"] = filename

            # Update the progress bar to reflect spectrum load in size (estimated by converting the dictionary to a string and checking its length)
            progress.update(len(str(spectrum)))

            # Yield the current spectrum, allowing the function to be used as a generator
            yield spectrum

        # Close the progress bar after all spectra have been loaded and yielded.
        progress.close()
