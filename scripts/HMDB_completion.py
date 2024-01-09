import concurrent.futures
from tqdm import tqdm
import pandas as pd
import os
import re

global HMDB_df
HMDB_df = pd.read_csv("../datas/HMDB.csv", sep=";", encoding="UTF-8")

def completing_spectrum(xml_path):
    """
    Modify an XML file by completing its spectrum with additional tags from a database.

    :param xml_path: Path to the XML file.
    :return: None if no modifications were made to the XML file, otherwise the modified XML file.
    """
    with open(xml_path, "r", encoding="UTF-8") as buffer:
        xml_content = buffer.read()

    HMDB_ID = re.search("<database-id>(.*)</database-id>", xml_content)  # hmdb id retrieval for csv matching
    if HMDB_ID:
        HMDB_ID = HMDB_ID.group(1)

        line = HMDB_df.loc[HMDB_df['EXTERNAL_ID'] == HMDB_ID]

        if not line.empty:
            INCHIKEY = line["INCHIKEY"].values[0]
            SMILES = line["SMILES"].values[0]
            INCHI = line["INCHI"].values[0]
            MOLECULAR_FORMULA = line["MOLECULAR_FORMULA"].values[0]
            ACC_MASS = line["ACC_MASS"].values[0]
        else:
            return None

        if re.search("<sample-mass>(.*)</sample-mass>", xml_content):
            old_tag = "</sample-mass>\n"
            new_tag = f"</sample-mass>\n  <inchikey>{INCHIKEY}</inchikey>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "</sample-mass>\n"
            new_tag = f"</sample-mass>\n  <SMILES>{SMILES}</SMILES>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "</sample-mass>\n"
            new_tag = f"</sample-mass>\n  <INCHI>{INCHI}</INCHI>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "</sample-mass>\n"
            new_tag = f"</sample-mass>\n  <MOLECULAR_FORMULA>{MOLECULAR_FORMULA}</MOLECULAR_FORMULA>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "</sample-mass>\n"
            new_tag = f"</sample-mass>\n  <ACC_MASS>{ACC_MASS}</ACC_MASS>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

        elif re.search("<sample-mass nil=\"true\"/>", xml_content):
            old_tag = "<sample-mass nil=\"true\"/>\n"
            new_tag = f"<sample-mass nil=\"true\"/>\n  <inchikey>{INCHIKEY}</inchikey>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "<sample-mass nil=\"true\"/>\n"
            new_tag = f"<sample-mass nil=\"true\"/>\n  <SMILES>{SMILES}</SMILES>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "<sample-mass nil=\"true\"/>\n"
            new_tag = f"<sample-mass nil=\"true\"/>\n  <INCHI>{INCHI}</INCHI>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "<sample-mass nil=\"true\"/>\n"
            new_tag = f"<sample-mass nil=\"true\"/>\n  <MOLECULAR_FORMULA>{MOLECULAR_FORMULA}</MOLECULAR_FORMULA>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

            old_tag = "<sample-mass nil=\"true\"/>\n"
            new_tag = f"<sample-mass nil=\"true\"/>\n  <ACC_MASS>{ACC_MASS}</ACC_MASS>\n"
            xml_content = xml_content.replace(old_tag, new_tag)

        with open(xml_path, "w", encoding="UTF-8") as buffer:
            buffer.write(xml_content)
    else:
        return None

def generate_xml_paths_list():
    """
    Generates a list of absolute file paths for all XML files in the specified directory.

    :return: A list of absolute file paths for all XML files.

    :rtype: list of str
    """
    # Liste pour stocker les chemins absolus
    xml_files_absolute_paths = []

    walk_dirs = os.walk("../INPUT/XML")  # Réinitialiser le générateur

    print("{:>80}".format("-- COMPLETING HMDB SPECTRUMS --"))

    for root, dirs, files in walk_dirs:
        for file in files:
            if file.endswith(".xml"):
                absolute_path = os.path.join(root, file)
                xml_files_absolute_paths.append(absolute_path)

    return xml_files_absolute_paths


def hmdb_completion_processing(xml_files_absolute_paths):
    """
    :param xml_files_absolute_paths: A list of absolute file paths to XML files containing HMDB spectra.
    :return: None

    This method processes the HMDB spectra by dividing the list of XML files into chunks and completing the spectra using multithreading. The spectra are completed by calling the `comple
    *ting_spectrum` function for each file in the chunk. The progress is tracked using a progress bar which displays the number of files completed.

    Note: The method does not return any value.
    """

    chunk_size = 5000
    progress_bar = tqdm(total=len(xml_files_absolute_paths), unit=" files", colour="green", desc="{:>80}".format("completing hmdb spectrums"))

    # Dividing the spectrum list into chunks
    for i in range(0, len(xml_files_absolute_paths), chunk_size):
        with concurrent.futures.ThreadPoolExecutor() as executor:
            chunk = xml_files_absolute_paths[i:i + chunk_size]
            executor.map(completing_spectrum, chunk)
            progress_bar.update(len(chunk))

    progress_bar.close()

xml_files_absolute_paths = generate_xml_paths_list()



