from rdkit.Chem.Descriptors import ExactMolWt
import concurrent.futures
import pubchempy as pcp
from rdkit import Chem
from tqdm import tqdm
import pandas as pd
from pycdk.pycdk import *
import uuid
import re
import os


def concatenate_clean_msp(clean_msp_path):
    CONCATENATE_LIST = []
    # append all spectrums of all cleaned files into CONCATENATE_LIST
    for files in os.listdir(clean_msp_path):
        if files.endswith("_clean.msp"):
            with open(os.path.join(clean_msp_path,files),"r",encoding="UTF-8") as buffer:
                temp = [element for element in buffer.read().split("\n\n") if element != "\n"]

            CONCATENATE_LIST.extend(temp)

    return CONCATENATE_LIST

def msp_to_csv(CONCATENATE_LIST):
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
    with open(msp_path, "r", encoding="UTF-8") as msp_buffer:
        content = msp_buffer.read()

    content = re.sub("charge: -\n","charge: -1\n",content,flags=re.I)

    with open(msp_path, "w", encoding="UTF-8") as msp_buffer:
        msp_buffer.write(content)

# def names_completion(CONCATENATE_LIST):
#     inchi_names = {}
#     inchi = "None"
#     name = "None"
#
#     for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
#         if re.search("INCHI: (.*)\n", spectrum):
#             inchi = re.search("INCHI: (.*)\n", spectrum).group(1)
#         if re.search("\nNAME: (.*)\n", spectrum):
#             name = re.search("\nNAME: (.*)\n", spectrum).group(1)
#
#         # If InChIKey and name are valid, add them to the dictionary
#         if inchi and name != "None":
#             inchi_names[inchi] = name
#
#     # Update missing names with corresponding names in dictionary list
#     updated_spetcrum_list = []
#     for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t    updating"):
#         if re.search("INCHI: (.*)\n", spectrum):
#             inchi = re.search("INCHI: (.*)\n", spectrum).group(1)
#         if re.search("\nNAME: (.*)\n", spectrum):
#             name = re.search("\nNAME: (.*)\n", spectrum).group(1)
#         if name == "None":
#             if inchi in inchi_names.keys():
#                 spectrum = re.sub("\nNAME: (.*)\n", f"\nNAME: {inchi_names[inchi]}\n", spectrum)
#
#         updated_spetcrum_list.append(spectrum)
#
#     return updated_spetcrum_list

# def names_completion(CONCATENATE_DF):
#     # Remplacer les 'None' par des chaînes vides dans la colonne "NAME"
#     CONCATENATE_DF['NAME'] = CONCATENATE_DF['NAME'].replace('None', '')
#
#     # Regrouper les lignes par "INCHI"
#     grouped = CONCATENATE_DF.groupby('INCHI')
#
#     # Créer une barre de progression tqdm
#     progress_bar = tqdm(CONCATENATE_DF, total=len(CONCATENATE_DF), unit=" rows", colour="green", desc="\t  processing")
#
#     # Itérer à travers les groupes
#     for _, group in grouped:
#         # Vérifier si le groupe a des noms vides
#         if group['NAME'].isnull().all():
#             # Obtenir le premier nom non nul du groupe
#             first_non_empty_name = group['NAME'].dropna().iloc[0]
#
#             # Si un nom non nul a été trouvé, attribuer ce nom à tout le groupe
#             if first_non_empty_name:
#                 group['NAME'] = first_non_empty_name
#
#         # Mettre à jour la barre de progression
#         progress_bar.update(1)
#
#     # Fermer la barre de progression
#     progress_bar.close()
#
#     # Mettre à jour le DataFrame original
#     CONCATENATE_DF.update(grouped.sum())
#
#     # Remplacer les valeurs vides par 'None' dans la colonne "NAME"
#     CONCATENATE_DF['NAME'].replace('', 'None', inplace=True)
#
#     return CONCATENATE_DF

def update_names(group):
    if group['NAME'].notnull().any():
        non_empty_name = group.loc[group['NAME'].notnull(), 'NAME'].iloc[0]
        group['NAME'] = non_empty_name
    return group

def names_completion(CONCATENATE_DF):
    CONCATENATE_DF['NAME'] = CONCATENATE_DF['NAME'].replace('None', '')

    grouped = CONCATENATE_DF.groupby('INCHI')
    updated_df = grouped.apply(update_names)

    # Remplacer les valeurs vides par 'None' dans la colonne "NAME"
    updated_df['NAME'] = updated_df['NAME'].replace('', 'None')

    return updated_df

def inchi_smiles_completion(CONCATENATE_LIST):
    inchikey_inchi = {}
    inchikey_smiles = {}

    inchikey = "None"
    inchi = "None"
    smiles = "None"

    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
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
    for spectrum in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t    updating"):
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
    CONCATENATE_LIST_temp = []
    for spectrums in tqdm(CONCATENATE_LIST, total=len(CONCATENATE_LIST), unit=" spectrums", colour="green", desc="\t  processing"):
        if not (re.search("SMILES: None\n",spectrums) or re.search("INCHI: None\n",spectrums)):
            CONCATENATE_LIST_temp.append(spectrums)

    return CONCATENATE_LIST_temp

def unique_id_generator():
    INPUT_path = "../INPUT/MSP"
    for files in os.listdir(INPUT_path):
        if files.endswith(".msp"):
            print(files)
            temp_list = []

            with open(os.path.join(INPUT_path,files),"r",encoding="UTF-8") as buffer:
                content = buffer.read()

            if not re.search("FRAGHUBID: (.*)\n",content):
                content = content.split("\n\n")

                for spectrums in tqdm(content, total=len(content), unit=" spectrums", colour="green", desc="\t  processing"):
                    spectrums = "FRAGHUBID: "+str(uuid.uuid4())+"\n"+spectrums
                    spectrums = re.sub("\n{2,}","\n",spectrums)
                    temp_list.append(spectrums)

                with open(os.path.join(INPUT_path,files),"w",encoding="UTF-8") as buffer:
                    buffer.write("\n\n\n".join(temp_list))

def apply_transformations(row):
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
    total_rows = len(CONCATENATE_DF)
    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="\t  generating")

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
    total_rows = len(CONCATENATE_DF)
    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="\t  processing")

    CONCATENATE_DF = CONCATENATE_DF.apply(mass_calculation, axis=1)
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return CONCATENATE_DF







