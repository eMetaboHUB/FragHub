from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from concurrent.futures import ThreadPoolExecutor
from rdkit import RDLogger, Chem
from tqdm import tqdm
import pandas as pd
import re
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

# Dossier contenant les fichiers CSV
folder_path = '../datas/pubchem_datas'
# Liste pour stocker chaque DataFrame
all_dfs = []
# Fonction pour lire un fichier CSV
def read_csv(file_path):
    return pd.read_csv(file_path, sep=';', quotechar='"', encoding='utf-8')

# Utiliser ThreadPoolExecutor pour lire les fichiers en parallèle
with ThreadPoolExecutor() as executor:
    futures = []
    for file_name in os.listdir(folder_path):
        if file_name.endswith('.csv'):
            file_path = os.path.join(folder_path, file_name)
            futures.append(executor.submit(read_csv, file_path))
    # Récupérer les résultats des futures
    for future in futures:
        all_dfs.append(future.result())

# Concaténer tous les DataFrames
global pubchem_datas
pubchem_datas = pd.concat(all_dfs, ignore_index=True)

inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey
indigo_smiles_correction_pattern = re.compile(r"\|[\s\S]*")

def apply_transformations(inchi_smiles):
    """
    Apply transformations to a given InChI or SMILES string.
    :param inchi_smiles: The InChI or SMILES string.
    :return: A dictionary containing the transformed values of the input string.
    """
    # Initialize an empty dictionary for transformed values
    transforms = {}

    # Check if the input string does not contain 'InChI='. If true, some corrections will be applied to the string
    if 'InChI=' not in inchi_smiles:
        inchi_smiles = re.sub(indigo_smiles_correction_pattern, "", inchi_smiles)

    # Check if the modified input string is indeed a string
    if isinstance(inchi_smiles, str):
        # Depending on whether it is in InChI format or SMILES format, appropriate conversion will be applied
        mol = Chem.MolFromInchi(inchi_smiles) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(inchi_smiles)

        # Checking if the conversion produced a valid result. If so, chemical information is extracted
        if mol is not None:
            transforms = {
                'INCHI': Chem.MolToInchi(mol),
                'INCHIKEY': Chem.MolToInchiKey(mol),
                'SMILES': Chem.MolToSmiles(mol),
                'FORMULA': CalcMolFormula(mol),
            }
        else:
            transforms = {
                'INCHI': '',
                'INCHIKEY': '',
                'SMILES': '',
                'FORMULA': '',
            }

        # Calculating the masses. If any errors occur during this process, default to blank strings
        if transforms:
            mol = Chem.MolFromInchi(transforms['INCHI']) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(transforms['SMILES'])
            if mol is not None:
                try:
                    transforms['EXACTMASS'] = ExactMolWt(mol)
                    transforms['AVERAGEMASS'] = MolWt(mol)
                except:
                    transforms['EXACTMASS'] = ''
                    transforms['AVERAGEMASS'] = ''
                    return transforms
            else:
                transforms['EXACTMASS'] = ''
                transforms['AVERAGEMASS'] = ''

    # Return the transformations dictionary
    return transforms

def map_transformations(row, unique_transforms):
    """
    Transforms the given row based on the unique transformations specified.
    :param row: A dictionary representing a row of data.
    :param unique_transforms: A dictionary of unique transformations.
    :return: The transformed row.
    """
    # Set original_inchi to the value of 'INCHI' in the row or None if it isn't available.
    original_inchi = row['INCHI'] if pd.notna(row['INCHI']) else None

    # Set original_smiles to the value of 'SMILES' in the row or None if it isn't available.
    original_smiles = row['SMILES'] if pd.notna(row['SMILES']) else None

    # If original_inchi is present and also there in unique_transforms, apply the transformations in the row.
    if original_inchi and original_inchi in unique_transforms:
        for key, value in unique_transforms[original_inchi].items():
            row[key] = value

    # If original_smiles is present and also there in unique_transforms, apply the transformations in the row.
    elif original_smiles and original_smiles in unique_transforms:
        for key, value in unique_transforms[original_smiles].items():
            row[key] = value

    # Return the transformed row
    return row

def merge_datas_from_pubchem(CONCATENATE_DF):
    # Premier merge par INCHIKEY
    CONCATENATE_DF = pd.merge(CONCATENATE_DF, pubchem_datas, on='INCHIKEY', how='left', suffixes=('', '_pubchem'))

    # Supprimer les colonnes suffixées après l'assignation
    CONCATENATE_DF.drop(columns=[col for col in CONCATENATE_DF.columns if '_pubchem' in col], inplace=True)

    # Masque pour les lignes où 'calculation' n'est pas 'done'
    not_done_mask = CONCATENATE_DF['calculation'] != 'done'

    # Deuxième fusion par INCHI sur les lignes non encore marquées comme 'done'
    CONCATENATE_DF.loc[not_done_mask] = pd.merge(
        CONCATENATE_DF[not_done_mask],
        pubchem_datas,
        on='INCHI',
        how='left',
        suffixes=('', '_pubchem')
    )

    # Supprimer à nouveau les colonnes suffixées après l'assignation
    CONCATENATE_DF.drop(columns=[col for col in CONCATENATE_DF.columns if '_pubchem' in col], inplace=True)

    # Mettre à jour le masque pour les lignes restantes
    not_done_mask = CONCATENATE_DF['calculation'] != 'done'

    # Troisième fusion par SMILES sur les lignes non encore marquées comme 'done'
    CONCATENATE_DF.loc[not_done_mask] = pd.merge(
        CONCATENATE_DF[not_done_mask],
        pubchem_datas,
        on='SMILES',
        how='left',
        suffixes=('', '_pubchem')
    )

    # Supprimer les colonnes suffixées après la troisième fusion
    CONCATENATE_DF.drop(columns=[col for col in CONCATENATE_DF.columns if '_pubchem' in col], inplace=True)

    return CONCATENATE_DF

def mols_derivation_and_calculation(CONCATENATE_DF):
    """
    Derives and calculates molecular properties based on unique INCHI and SMILES in the given DataFrame.
    :param CONCATENATE_DF: DataFrame containing the INCHI and SMILES columns.
    :return: DataFrame with calculated molecular properties.
    """
    # Initialize the calculation column
    CONCATENATE_DF['calculation'] = ''

    # Merge data from pubchem
    CONCATENATE_DF = merge_datas_from_pubchem(CONCATENATE_DF)

    # Filter for rows where calculation is not done
    not_done_mask = CONCATENATE_DF['calculation'] != 'done'

    # Concatenate 'INCHI' and 'SMILES' columns, drop all null entries and get unique values from result
    unique_inchi_smiles = pd.concat([
        CONCATENATE_DF.loc[not_done_mask, 'INCHI'],
        CONCATENATE_DF.loc[not_done_mask, 'SMILES']
    ]).dropna().unique()

    # For each unique inchi_smiles run 'apply_transformations' and store result to a dictionary
    unique_transforms = {inchi_smiles: apply_transformations(inchi_smiles) for inchi_smiles in tqdm(
        unique_inchi_smiles, unit=" rows", colour="green", desc="{:>70}".format("derivation and calculation"))}

    # Initialize the tqdm progress bar for the coming dataframe calculations
    tqdm.pandas(unit="rows", colour="green", desc="{:>70}".format("updating dataframe"))

    # Apply the transformations stored in 'unique_transforms' to the dataframe 'CONCATENATE_DF'
    CONCATENATE_DF = CONCATENATE_DF.progress_apply(map_transformations, axis=1, args=(unique_transforms,))

    # Check if 'INCHIKEY' in 'CONCATENATE_DF' dataframe matches 'inchikey_pattern' string pattern
    mask = CONCATENATE_DF['INCHIKEY'].str.fullmatch(inchikey_pattern, na=False)

    # Apply the mask to the dataframe
    CONCATENATE_DF = CONCATENATE_DF[mask]

    # Drop null entries in 'EXACTMASS', 'AVERAGEMASS', 'SMILES', 'INCHI', and 'INCHIKEY' columns of 'CONCATENATE_DF'
    CONCATENATE_DF = CONCATENATE_DF.dropna(subset=['EXACTMASS', 'AVERAGEMASS', 'SMILES', 'INCHI', 'INCHIKEY'])

    # Return the transformed DataFrame
    return CONCATENATE_DF
