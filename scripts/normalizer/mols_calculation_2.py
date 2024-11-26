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
            mol = Chem.MolFromInchi(transforms['INCHI']) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(
                transforms['SMILES'])
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


def complete_from_pubchem_datas(CONCATENATE_DF, pubchem_datas):
    """
    Complète le DataFrame CONCATENATE_DF avec des informations provenant de pubchem_datas en utilisant plusieurs fusions successives.

    :param CONCATENATE_DF: DataFrame contenant les colonnes INCHIKEY, INCHI, SMILES, FORMULA, NAME, EXACTMASS, AVERAGEMASS.
    :param pubchem_datas: DataFrame contenant des informations supplémentaires.
    :return: Tuple contenant deux DataFrames :
             - DataFrame CONCATENATE_DF complété avec `done`.
             - DataFrame CONCATENATE_DF non complété.
    """
    # Ajoute une colonne pour suivre l'état de calcul
    CONCATENATE_DF['calculation'] = ''

    # Définit une fonction pour la fusion et mise à jour du statut 'done'
    def merge_and_mark_done(df1, df2, key_col):
        merged_df = pd.merge(df1, df2, on=key_col, how='left', suffixes=('', '_pubchem'))
        for col in ['INCHI', 'SMILES', 'FORMULA', 'NAME', 'EXACTMASS', 'AVERAGEMASS']:
            if col + '_pubchem' in merged_df.columns:
                # Complète les colonnes manquantes
                merged_df[col].fillna(merged_df[col + '_pubchem'], inplace=True)
                merged_df.drop(columns=[col + '_pubchem'], inplace=True)

        # Marque les lignes fusionnées avec 'done'
        merged_df.loc[merged_df[key_col].notna(), 'calculation'] = 'done'
        return merged_df

    # Premier merge sur 'INCHIKEY'
    CONCATENATE_DF = merge_and_mark_done(CONCATENATE_DF, pubchem_datas, 'INCHIKEY')

    # Deuxième merge sur 'INCHI' pour les lignes qui ne sont pas encore 'done'
    not_done_df = CONCATENATE_DF[CONCATENATE_DF['calculation'] != 'done']
    updated_df = merge_and_mark_done(not_done_df, pubchem_datas, 'INCHI')
    CONCATENATE_DF.update(updated_df)

    # Troisième merge sur 'SMILES' pour les lignes qui ne sont pas encore 'done'
    not_done_df = CONCATENATE_DF[CONCATENATE_DF['calculation'] != 'done']
    updated_df = merge_and_mark_done(not_done_df, pubchem_datas, 'SMILES')
    CONCATENATE_DF.update(updated_df)

    # Séparation des DataFrames basée sur 'calculation' état
    done_df = CONCATENATE_DF[CONCATENATE_DF['calculation'] == 'done'].copy().drop(columns=['calculation'])
    not_done_df = CONCATENATE_DF[CONCATENATE_DF['calculation'] != 'done'].copy().drop(columns=['calculation'])

    return done_df, not_done_df


def mols_derivation_and_calculation(CONCATENATE_DF):
    """
    Derives and calculates molecular properties based on unique INCHI and SMILES in the given DataFrame.
    :param CONCATENATE_DF: DataFrame containing the INCHI and SMILES columns.
    :return: DataFrame with calculated molecular properties.
    """
    done_df, not_done_df = complete_from_pubchem_datas(CONCATENATE_DF, pubchem_datas)

    del CONCATENATE_DF

    # Concaténe les colonnes 'INCHI' et 'SMILES', élimine les entrées nulles et obtient des valeurs uniques
    unique_inchi_smiles = pd.concat([not_done_df['INCHI'], not_done_df['SMILES']]).dropna().unique()

    # Applique 'apply_transformations' à chaque valeur unique de inchi_smiles et stocke le résultat dans un dictionnaire
    unique_transforms = {inchi_smiles: apply_transformations(inchi_smiles) for inchi_smiles in
                         tqdm(unique_inchi_smiles, unit=" rows", colour="green",
                              desc="{:>70}".format("derivation and calculation"))}

    # Initialize the tqdm progress bar for the coming dataframe calculations
    tqdm.pandas(unit="rows", colour="green", desc="{:>70}".format("updating dataframe"))

    # Apply the transformations stored in 'unique_transforms' to the dataframe 'not_done_df'
    not_done_df = not_done_df.progress_apply(map_transformations, axis=1, args=(unique_transforms,))

    # Check if 'INCHIKEY' in 'not_done_df' dataframe matches 'inchikey_pattern' string pattern
    mask = not_done_df['INCHIKEY'].str.fullmatch(inchikey_pattern, na=False)

    # Apply the mask to the dataframe
    not_done_df = not_done_df[mask]

    # Drop null entries in 'EXACTMASS', 'AVERAGEMASS', 'SMILES', 'INCHI', and 'INCHIKEY' columns of 'not_done_df'
    not_done_df = not_done_df.dropna(subset=['EXACTMASS', 'AVERAGEMASS', 'SMILES', 'INCHI', 'INCHIKEY'])

    # Concaténer les dataframes done et not_done après transformation
    CONCATENATE_DF = pd.concat([done_df, not_done_df], ignore_index=True)

    del done_df
    del not_done_df

    # Return the transformed DataFrame
    return CONCATENATE_DF
