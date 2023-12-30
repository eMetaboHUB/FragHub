from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit import RDLogger, Chem
from pycdk.pycdk import *
from tqdm import tqdm
import pandas as pd

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
        if transforms:
            mol = MolFromInchi(transforms['INCHI']) if 'InChI=' in inchi_smiles else MolFromSmiles(transforms['SMILES'])
            if mol is not None:
                try:
                    transforms['EXACTMASS'] = str(getMolExactMass(mol))
                    transforms['AVERAGEMASS'] = str(getMolNaturalMass(mol))
                except:
                    return transforms

    return transforms


def map_transformations(row, unique_transforms):
    """
    Maps transformations to each row based on the values of 'INCHI' or 'SMILES' columns.

    :param row: A pandas DataFrame row.
    :param unique_transforms: A dictionary that maps unique values of 'INCHI' or 'SMILES' to transformations.
    :return: The updated row with the transformations applied.
    """
    if row['INCHI'] in unique_transforms:
        row.update(pd.Series(unique_transforms[row['INCHI']]))
    elif row['SMILES'] in unique_transforms:
        row.update(pd.Series(unique_transforms[row['SMILES']]))
    return row


def mols_derivation_and_calculation(CONCATENATE_DF):
    """
    Derives and calculates molecular properties based on unique INCHI and SMILES in the given DataFrame.

    :param CONCATENATE_DF: DataFrame containing the INCHI and SMILES columns.
    :return: DataFrame with calculated molecular properties.
    """
    unique_inchi_smiles = pd.concat([CONCATENATE_DF['INCHI'], CONCATENATE_DF['SMILES']]).dropna().unique()

    # Creating a dict that maps each unique INCHI or SMILES to its respective transformations.
    unique_transforms = {inchi_smiles: apply_transformations(inchi_smiles) for inchi_smiles in tqdm(unique_inchi_smiles, unit=" rows", colour="green", desc="\t  generating")}

    # Using apply to apply the transformations
    CONCATENATE_DF = CONCATENATE_DF.apply(map_transformations, axis=1, args=(unique_transforms,))

    return CONCATENATE_DF