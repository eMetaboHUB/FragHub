from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger
from pycdk.pycdk import *
from tqdm import tqdm
import pandas as pd
from rdkit import Chem

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

def apply_transformations(row):
    unique_transforms = {}

    if not pd.isna(row['INCHI']):
        mol = Chem.MolFromInchi(row['INCHI'])
        if mol is not None:
            unique_transforms[row['INCHI']] = {
                'INCHI': Chem.MolToInchi(mol),
                'INCHIKEY': Chem.MolToInchiKey(mol),
                'SMILES': Chem.MolToSmiles(mol),
                'FORMULA': CalcMolFormula(mol),
            }
    elif not pd.isna(row['SMILES']):
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is not None:
            unique_transforms[row['SMILES']] = {
                'INCHI': Chem.MolToInchi(mol),
                'INCHIKEY': Chem.MolToInchiKey(mol),
                'SMILES': Chem.MolToSmiles(mol),
                'FORMULA': CalcMolFormula(mol),
            }
    return unique_transforms

def mols_derivator(CONCATENATE_DF):
    """
    :param CONCATENATE_DF: DataFrame containing the data to be transformed.
    :return: DataFrame with the transformed data.

    This method applies transformations to each row of the input DataFrame. It tracks the progress using a progress bar provided by the `tqdm` library.
    The transformations are applied by calling the `apply_transformations` function on each row. After all transformations have been collected, they are applied to DataFrame. The updated DataFrame is returned.
    """

    unique_transforms = {}
    total_rows = len(CONCATENATE_DF)

    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="\t  generating")

    for i in range(total_rows):
        row = CONCATENATE_DF.iloc[i]
        unique_transforms.update(apply_transformations(row))
        t.update()

    t.close()

    for key, val in unique_transforms.items():
        CONCATENATE_DF.loc[CONCATENATE_DF['INCHI'] == key] = val
        CONCATENATE_DF.loc[CONCATENATE_DF['SMILES'] == key] = val

    return CONCATENATE_DF

def mass_calculation(row):
    """
    Mass calculation method.

    :param row: pandas DataFrame row containing 'SMILES' column
    :return: pandas DataFrame row with calculated 'EXACTMASS' and 'AVERAGEMASS' columns

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
    Calculate mass for each row in the given DataFrame.

    :param CONCATENATE_DF: The DataFrame containing the data.
    :return: The DataFrame with mass values calculated and updated.
    """
    total_rows = len(CONCATENATE_DF)
    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="\t  processing")

    CONCATENATE_DF = CONCATENATE_DF.apply(mass_calculation, axis=1)
    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

    return CONCATENATE_DF