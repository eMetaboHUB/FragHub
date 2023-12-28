from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger
from pycdk.pycdk import *
from tqdm import tqdm
import pandas as pd
from rdkit import Chem

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

def apply_transformations(row):
    """
    Apply transformations to the given row.

    :param row: The row to apply transformations to.
    :return: The row with transformed values.

    """
    if not pd.isna(row['INCHI']):
        mol = Chem.MolFromInchi(row['INCHI'])
        if mol is not None:
            row['INCHI'] = Chem.MolToInchi(mol)
            row['INCHIKEY'] = Chem.MolToInchiKey(mol)
            row['SMILES'] = Chem.MolToSmiles(mol)
            row['FORMULA'] = CalcMolFormula(mol)
    elif not pd.isna(row['SMILES']):
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is not None:
            row['SMILES'] = Chem.MolToSmiles(mol)
            row['INCHI'] = Chem.MolToInchi(mol)
            row['INCHIKEY'] = Chem.MolToInchiKey(mol)
            row['FORMULA'] = CalcMolFormula(mol)
    return row

def mols_derivator(CONCATENATE_DF):
    """
    :param CONCATENATE_DF: DataFrame containing the data to be transformed.
    :return: DataFrame with the transformed data.

    This method applies transformations to each row of the input DataFrame in parallel using the `apply` function. It tracks the progress using a progress bar provided by the `tqdm` library
    *. The transformations are applied by calling the `apply_transformations` function on each row. After all transformations are applied, the updated DataFrame is returned.

    Example usage:
    ```
    result = mols_derivator(input_df)
    ```
    """
    total_rows = len(CONCATENATE_DF)
    t = tqdm(total=total_rows, unit=" rows", colour="green", desc="\t  generating")

    CONCATENATE_DF = CONCATENATE_DF.apply(apply_transformations, axis=1)

    t.update(total_rows)

    # Fermer la barre de progression
    t.close()

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