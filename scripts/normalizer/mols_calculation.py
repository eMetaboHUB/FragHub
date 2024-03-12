from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit import RDLogger, Chem
from tqdm import tqdm
import pandas as pd
import re

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey
indigo_smiles_correction_pattern = re.compile(r"\|[\s\S]*")

def apply_transformations(inchi_smiles):
    """
    Apply transformations to a given InChI or SMILES string.

    :param inchi_smiles: The InChI or SMILES string.
    :return: A dictionary containing the transformed values of the input string.
    """
    transforms = {}

    if 'InChI=' not in inchi_smiles:
        inchi_smiles = re.sub(indigo_smiles_correction_pattern, "", inchi_smiles)

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
        else:
            transforms = {
                'INCHI': '',
                'INCHIKEY': '',
                'SMILES': '',
                'FORMULA': '',
            }
        # Mass calculation
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

    return transforms

def map_transformations(row, unique_transforms):
    """
    Transforms the given row based on the unique transformations specified.

    :param row: A dictionary representing a row of data.
    :param unique_transforms: A dictionary of unique transformations.
    :return: The transformed row.
    """
    original_inchi = row['INCHI'] if pd.notna(row['INCHI']) else None
    original_smiles = row['SMILES'] if pd.notna(row['SMILES']) else None

    if original_inchi and original_inchi in unique_transforms:
        for key, value in unique_transforms[original_inchi].items():
            row[key] = value
    elif original_smiles and original_smiles in unique_transforms:
        for key, value in unique_transforms[original_smiles].items():
            row[key] = value

    return row

def mols_derivation_and_calculation(CONCATENATE_DF):
    """
    Derives and calculates molecular properties based on unique INCHI and SMILES in the given DataFrame.

    :param CONCATENATE_DF: DataFrame containing the INCHI and SMILES columns.
    :return: DataFrame with calculated molecular properties.
    """
    unique_inchi_smiles = pd.concat([CONCATENATE_DF['INCHI'], CONCATENATE_DF['SMILES']]).dropna().unique()

    # Creating a dict that maps each unique INCHI or SMILES to its respective transformations.
    unique_transforms = {inchi_smiles: apply_transformations(inchi_smiles) for inchi_smiles in tqdm(unique_inchi_smiles, unit=" rows", colour="green", desc="{:>70}".format("derivation and calculation"))}

    # Set up progress apply with tqdm
    tqdm.pandas(unit=" rows", colour="green", desc="{:>70}".format("updating dataframe"))

    # Using apply to apply the transformations
    CONCATENATE_DF = CONCATENATE_DF.progress_apply(map_transformations, axis=1, args=(unique_transforms,))

    mask = CONCATENATE_DF['INCHIKEY'].str.fullmatch(inchikey_pattern, na=False)
    CONCATENATE_DF = CONCATENATE_DF[mask]

    CONCATENATE_DF = CONCATENATE_DF.dropna(subset=['EXACTMASS', 'AVERAGEMASS', 'SMILES', 'INCHI', 'INCHIKEY'])

    return CONCATENATE_DF