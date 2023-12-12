import re


def load_spectrum_list(msp_file_path):
    """
    Load spectrum list from the MSP file.

    :param msp_file_path: The file path of the MSP file.
    :return: A list of spectrum strings.
    """
    with open(msp_file_path, 'r') as file:
        content = file.read()
    spectrum_list = content.split("\n\n")

    return spectrum_list

def extract_metadata_and_peak_list(spectrum):
    """
    Extracts metadata and peak list from a spectrum.

    :param spectrum: The spectrum containing metadata and peak list.
    :return: A tuple containing the extracted metadata and peak list.
    """
    if re.search("([\s\S]*:.[0-9]*\n)([(\d+\.?\d*)\s+(\d+\.?\d*)]*)",spectrum):
        match = re.search("([\s\S]*:.[0-9]*\n)([(\d+\.?\d*)\s+(\d+\.?\d*)]*)", spectrum)
        metadata, peak_list = match.group(1), match.group(2)

        return metadata, peak_list

def parse_metadata_and_peak_list(spectrum):
    metadata, peak_list = extract_metadata_and_peak_list(spectrum)




def parse_msp_file(msp_file_path):
    spectrum_list = load_spectrum_list(msp_file_path)

    for spectrum in spectrum_list:
        metadata,peak_list = parse_metadata_and_peak_list(spectrum)

