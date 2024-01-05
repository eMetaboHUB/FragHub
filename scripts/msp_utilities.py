from rdkit.Chem.Descriptors import ExactMolWt
import concurrent.futures
from pycdk.pycdk import *
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
import uuid
import re
import os


def concatenate_clean_msp(clean_msp_path):
    """
    Concatenates all spectrums of all cleaned files in the given directory.

    :param clean_msp_path: Path to the directory containing the cleaned files.
    :return: A list of concatenated spectrums.
    """
    CONCATENATE_LIST = []
    # append all spectrums of all cleaned files into CONCATENATE_LIST
    for files in os.listdir(clean_msp_path):
        if files.endswith("_clean.msp"):
            with open(os.path.join(clean_msp_path,files),"r",encoding="UTF-8") as buffer:
                temp = [element for element in buffer.read().split("\n\n") if element != "\n"]

            CONCATENATE_LIST.extend(temp)

    return CONCATENATE_LIST

def msp_to_csv(CONCATENATE_LIST):
    """
    This method converts a list of spectra in MSP format to a CSV format.

    :param CONCATENATE_LIST: List of spectra in MSP format
    :return: Converted data in DataFrame format

    """
    dictionary = {}
    CONCATENATE_DF = pd.DataFrame()
    first = True
    empty = False
    if len(CONCATENATE_LIST) == 0:
        empty = True

    for spectrum in CONCATENATE_LIST:
        if spectrum != "\n":
            fields = re.findall(r"(.+?):(.*)\n", spectrum)
            if first:
                for element in fields:
                    dictionary[element[0]] = []
                dictionary["PEAKS_LIST"] = []
                first = False
            if not first:
                for element in fields:
                    dictionary[element[0]].append(element[1].strip(" "))

            if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2))

    if not empty:
        # Creating Dataframe
        CONCATENATE_DF = pd.DataFrame.from_dict(dictionary)
        CONCATENATE_DF = CONCATENATE_DF.drop("PARENTMASS", axis=1)
        CONCATENATE_DF.insert(12, "EXACTMASS", ["None" for i in range(len(dictionary["PEAKS_LIST"]))])
        CONCATENATE_DF.insert(13, "AVERAGEMASS", ["None" for i in range(len(dictionary["PEAKS_LIST"]))])

    return CONCATENATE_DF

def correct_uncomplete_charge(msp_path):
    """
    Corrects an incomplete charge in the MSP file.

    :param msp_path: The path to the MSP file.
    :return: None
    """
    with open(msp_path, "r", encoding="UTF-8") as msp_buffer:
        content = msp_buffer.read()

    content = re.sub("charge: -\n","charge: -1\n",content,flags=re.I)

    with open(msp_path, "w", encoding="UTF-8") as msp_buffer:
        msp_buffer.write(content)

def names_completion(CONCATENATE_DF):
    """
    Perform name completion for a given DataFrame.

    :param CONCATENATE_DF: The DataFrame containing the 'NAME' column to perform name completion on.
    :return: The modified DataFrame with name completion applied.

    """
    # Définir le nom de la barre de progression
    tqdm.pandas(total=len(CONCATENATE_DF), colour="green", unit=" row", desc="{:>40}".format("updating names"))

    # Appliquer la transformation par groupe avec une barre de progression
    CONCATENATE_DF['NAME'] = CONCATENATE_DF.groupby('INCHI')['NAME'].progress_transform(lambda group: group.fillna(group.dropna().iloc[0] if group.dropna().size > 0 else ''))

    return CONCATENATE_DF

def inchi_smiles_completion(CONCATENATE_LIST):
    """
    This method updates missing inchi and smiles values in a list of spectra using corresponding inchikey values.

    :param CONCATENATE_LIST: A list of spectra to be updated
    :return: The updated list of spectra
    """
    inchikey_inchi = {}
    inchikey_smiles = {}

    inchikey = "None"
    inchi = "None"
    smiles = "None"

    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="{:>40}".format("processing")):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        if re.search("\nINCHI: (.*)\n", spectrum):
            inchi = re.search("\nINCHI: (.*)\n", spectrum).group(1)
        if re.search("\nSMILES: (.*)\n", spectrum):
            smiles = re.search("\nSMILES: (.*)\n", spectrum).group(1)

        if inchikey and inchi != "None":
            inchikey_inchi[inchikey] = inchi

        if inchikey and smiles != "None":
            inchikey_smiles[inchikey] = smiles

    # Update missing inchi/smiles with corresponding inchikey in dictionary list
    updated_spetcrum_list = []
    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="{:>40}".format("updating")):
        if re.search("INCHIKEY: (.*)\n", spectrum):
            inchikey = re.search("INCHIKEY: (.*)\n", spectrum).group(1)
        # INCHI
        if re.search("\nINCHI: (.*)\n", spectrum):
            inchi = re.search("\nINCHI: (.*)\n", spectrum).group(1)
        if inchi == "None":
            if inchikey in inchikey_inchi.keys():
                spectrum = re.sub("\nINCHI: (.*)\n", rf"\nINCHI: {re.escape(inchikey_inchi[inchikey])}\n", spectrum)
                spectrum = re.sub(r"\\\\",r"\\",spectrum)
        # SMILES
        if re.search("\nSMILES: (.*)\n", spectrum):
            smiles = re.search("\nSMILES: (.*)\n", spectrum).group(1)
        if smiles == "None":
            if inchikey in inchikey_smiles.keys():
                spectrum = re.sub("\nSMILES: (.*)\n", rf"\nSMILES: {re.escape(inchikey_smiles[inchikey])}\n", spectrum)
                spectrum = re.sub(r"\\\\", r"\\", spectrum)

        updated_spetcrum_list.append(spectrum)

    return updated_spetcrum_list

def remove_no_smiles_inchi(CONCATENATE_LIST):
    """
    :param CONCATENATE_LIST: A list of strings representing spectrums.
    :return: A new list containing only spectrums that have both "SMILES: None" and "INCHI: None" strings removed.
    """
    CONCATENATE_LIST_temp = []
    for spectrums in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="{:>40}".format("processing")):
        if not (re.search("SMILES: None\n",spectrums) or re.search("INCHI: None\n",spectrums)):
            CONCATENATE_LIST_temp.append(spectrums)

    return CONCATENATE_LIST_temp

def unique_id_generator():
    """
    Generates unique IDs for each spectrum file found in the INPUT_path directory.

    :return: None
    """
    INPUT_path = "../INPUT/MSP"
    for files in os.listdir(INPUT_path):
        if files.endswith(".msp"):
            print(files)
            temp_list = []

            with open(os.path.join(INPUT_path,files),"r",encoding="UTF-8") as buffer:
                content = buffer.read()

            if not re.search("FRAGHUBID: (.*)\n",content):
                content = content.split("\n\n")
                content = [spectrum for spectrum in content if len(spectrum) > 10]

                for spectrums in tqdm(content, total=len(content), unit=" spectrums", colour="green", desc="{:>40}".format("processing")):
                    spectrums = "FRAGHUBID: "+str(uuid.uuid4())+"\n"+spectrums
                    spectrums = re.sub("\n{2,}","\n",spectrums)
                    temp_list.append(spectrums)

                with open(os.path.join(INPUT_path,files),"w",encoding="UTF-8") as buffer:
                    buffer.write("\n\n\n".join(temp_list))

def apply_transformations(row):
    """
    Apply transformations to a given row.

    :param row: A row containing chemical information.
    :return: The modified row with updated chemical information.

    The `apply_transformations` method takes a row as input and applies chemical transformations to the row's INCHI and SMILES values. If the INCHI value is not NaN, the method converts
    * it to a molecule object using the `MolFromInchi` function from the RDKit library. If a valid molecule object is obtained, the method updates the INCHI, INCHIKEY, and SMILES values
    * of the row using the corresponding conversion functions from the RDKit library. If the INCHI value is NaN and the SMILES value is not NaN, the method performs a similar transformation
    * using the `MolFromSmiles` function. Finally, the modified row with updated chemical information is returned.

    Note: This method assumes the use of the RDKit library and the presence of the pandas library.

    Example usage:

    ```
    row = {'INCHI': 'InChI String', 'SMILES': 'SMILES String'}
    modified_row = apply_transformations(row)
    print(modified_row)
    ```
    """
    if not pd.isna(row['INCHI']):
        mol = Chem.MolFromInchi(row['INCHI'])
        if mol is not None:
            row['INCHI'] = Chem.MolToInchi(mol)
            row['INCHIKEY'] = Chem.MolToInchiKey(mol)
            row['SMILES'] = Chem.MolToSmiles(mol)
    elif not pd.isna(row['SMILES']):
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is not None:
            row['SMILES'] = Chem.MolToSmiles(mol)
            row['INCHI'] = Chem.MolToInchi(mol)
            row['INCHIKEY'] = Chem.MolToInchiKey(mol)
    return row

def mols_derivator(CONCATENATE_DF):
    """
    :param CONCATENATE_DF: pandas DataFrame containing the data to be processed
    :return: pandas DataFrame with the processed data

    This method takes a pandas DataFrame as input and performs the following operations:
    1. Replaces 'None' values in the 'INCHI' column with NaN.
    2. Replaces 'None' values in the 'SMILES' column with NaN.
    3. Replaces 'None' values in the 'INCHIKEY' column with NaN.
    4. Applies the apply_transformations function to each row of the DataFrame.
    5. Drops rows where 'INCHI' column is NaN.
    6. Drops rows where 'SMILES' column is NaN.
    7. Drops rows where 'INCHIKEY' column is NaN.
    8. Updates the progress bar with the total number of rows.
    9. Closes the progress bar.
    10. Returns the modified DataFrame.
    """
    total_rows = len(CONCATENATE_DF)
    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="{:>40}".format("generating"))

    # supprimer les spectres sans InChI || SMILES || InChIKey
    CONCATENATE_DF['INCHI'] = CONCATENATE_DF['INCHI'].replace('None', float('nan'))
    CONCATENATE_DF['SMILES'] = CONCATENATE_DF['SMILES'].replace('None', float('nan'))
    CONCATENATE_DF['INCHIKEY'] = CONCATENATE_DF['INCHIKEY'].replace('None', float('nan'))

    CONCATENATE_DF = CONCATENATE_DF.apply(apply_transformations, axis=1)

    CONCATENATE_DF = CONCATENATE_DF.dropna(subset=['INCHI'])
    CONCATENATE_DF = CONCATENATE_DF.dropna(subset=['SMILES'])
    CONCATENATE_DF = CONCATENATE_DF.dropna(subset=['INCHIKEY'])

    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return CONCATENATE_DF

def mass_calculation(row):
    """
    Calculate the exact mass and average mass of a molecular structure.

    :param row: A pandas DataFrame row containing a 'SMILES' column.
    :return: A modified row with 'EXACTMASS' and 'AVERAGEMASS' columns.

    """
    if not pd.isna(row['SMILES']):
        SMILES = str(row['SMILES'])
        if SMILES != "None":
            try:
                mol = MolFromSmiles(SMILES)
                if mol is not None:
                    row['EXACTMASS'] = str(getMolExactMass(mol))
            except:
                return row
        # CDK average mass
        SMILES = str(row['SMILES'])
        if SMILES != "None":
            try:
                mol = MolFromSmiles(SMILES)
                if mol is not None:
                    row['AVERAGEMASS'] = str(getMolNaturalMass(mol))
            except:
                return row

    return row

def mass_calculator(CONCATENATE_DF):
    """
    Calculate the mass for each row in the given DataFrame.

    :param CONCATENATE_DF: The DataFrame containing the rows to calculate the mass for.
    :return: The updated DataFrame with the calculated mass values.
    """
    total_rows = len(CONCATENATE_DF)
    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="{:>40}".format("processing"))

    CONCATENATE_DF = CONCATENATE_DF.apply(mass_calculation, axis=1)
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return CONCATENATE_DF







