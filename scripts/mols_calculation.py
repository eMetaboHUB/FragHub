from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger
from pycdk.pycdk import *
from tqdm import tqdm
import pandas as pd
from rdkit import Chem

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

def apply_transformations(inchi_smiles):
    """
    Apply transformations to a given InChI or SMILES string.

    :param inchi_smiles: The InChI or SMILES string.
    :return: A dictionary containing the transformed values of the input string.
    """
    transforms = {}

    if isinstance(inchi_smiles, str):
        mol = Chem.MolFromInchi(inchi_smiles) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(inchi_smiles)
        # Mol harmonization
        if mol is not None:
            transforms = {
                'INCHI': Chem.MolToInchi(mol),
                'INCHIKEY': Chem.MolToInchiKey(mol),
                'SMILES': Chem.MolToSmiles(mol),
                'FORMULA': CalcMolFormula(mol),
            }
        # Mass calculation
        mol = MolFromInchi(inchi_smiles) if 'InChI=' in inchi_smiles else MolFromSmiles(inchi_smiles)
        if mol is not None:
            try:
                transforms['EXACTMASS'] = str(getMolExactMass(mol))
                transforms['AVERAGEMASS'] = str(getMolNaturalMass(mol))
            except:
                return transforms

    return transforms

def map_transformations(row, unique_transforms):
    """
    Transforms a DataFrame row using a dictionary of unique transformations.

    :param row: A row from a DataFrame.
    :param unique_transforms: A dictionary containing unique transformations, with keys as INCHI or SMILES and values as transformed data.
    :return: A transformed row as a pandas Series, if the row's INCHI or SMILES exists in the unique_transforms dictionary. Otherwise, returns the original row.

    """
    if row['INCHI'] in unique_transforms:
        return pd.Series(unique_transforms[row['INCHI']])
    elif row['SMILES'] in unique_transforms:
        return pd.Series(unique_transforms[row['SMILES']])
    else:
        return row

def mols_derivation_and_calculation(CONCATENATE_DF):
    """
    :param CONCATENATE_DF: A pandas DataFrame containing columns 'INCHI' and 'SMILES'.
    :return: The modified CONCATENATE_DF DataFrame with additional columns containing the respective transformations for each unique INCHI or SMILES.

    This method derives and calculates transformations for each unique INCHI or SMILES in the given CONCATENATE_DF DataFrame. It first identifies the unique INCHI or SMILES and creates a
    * dictionary that maps each unique INCHI or SMILES to its respective transformations. Then, it applies the transformations to each row of the CONCATENATE_DF DataFrame using the map_transform
    *ations function. Finally, it returns the modified CONCATENATE_DF DataFrame.
    """
    unique_inchi_smiles = pd.concat([CONCATENATE_DF['INCHI'], CONCATENATE_DF['SMILES']]).dropna().unique()

    # Creating a dict that maps each unique INCHI or SMILES to its respective transformations.
    unique_transforms = {inchi_smiles: apply_transformations(inchi_smiles) for inchi_smiles in tqdm(unique_inchi_smiles, unit=" rows", colour="green", desc="generating")}

    CONCATENATE_DF = CONCATENATE_DF.apply(map_transformations, axis=1, args=(unique_transforms,))

    return CONCATENATE_DF

# def mass_calculation(row):
#     """
#     Mass calculation method.
#
#     :param row: pandas DataFrame row containing 'SMILES' column
#     :return: pandas DataFrame row with calculated 'EXACTMASS' and 'AVERAGEMASS' columns
#
#     """
#     if not pd.isna(row['SMILES']):
#         SMILES = str(row['SMILES'])
#         if SMILES != "None":
#             try:
#                 mol = MolFromSmiles(SMILES)
#                 if mol is not None:
#                     row['EXACTMASS'] = str(getMolExactMass(mol))
#             except:
#                 return row
#         # CDK average mass
#         SMILES = str(row['SMILES'])
#         if SMILES != "None":
#             try:
#                 mol = MolFromSmiles(SMILES)
#                 if mol is not None:
#                     row['AVERAGEMASS'] = str(getMolNaturalMass(mol))
#             except:
#                 return row
#
#     return row
#
# def mass_calculator(CONCATENATE_DF):
#     """
#     Calculate mass for each row in the given DataFrame.
#
#     :param CONCATENATE_DF: The DataFrame containing the data.
#     :return: The DataFrame with mass values calculated and updated.
#     """
#     total_rows = len(CONCATENATE_DF)
#     t = tqdm(total=total_rows, unit=" rows", colour="green", desc="\t  processing")
#
#     CONCATENATE_DF = CONCATENATE_DF.apply(mass_calculation, axis=1)
#     t.update(total_rows)
#
#     # Fermer la barre de progression
#     t.close()
#
#     return CONCATENATE_DF