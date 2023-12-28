from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from pycdk.pycdk import *
from rdkit import Chem

def mols_derivator(metadata_dict):
    """
    Derives molecular properties from the given metadata dictionary.

    :param metadata_dict: A dictionary containing metadata of a molecule.
                          It should have the keys "INCHI" and "SMILES".
                          The values corresponding to these keys should be valid InChI and SMILES representations, respectively.
                          Note: At least one of the keys should have a value.

    :return: The metadata dictionary with additional derived properties.
             The keys "INCHI", "INCHIKEY", "SMILES", and "FORMULA" will be added with their respective derived values.

    """
    inchi = metadata_dict["INCHI"]
    smiles = metadata_dict["SMILES"]

    if inchi:
        mol = Chem.MolFromInchi(inchi)
        if mol is not None:
            metadata_dict['INCHI'] = Chem.MolToInchi(mol)
            metadata_dict['INCHIKEY'] = Chem.MolToInchiKey(mol)
            metadata_dict['SMILES'] = Chem.MolToSmiles(mol)
            metadata_dict['FORMULA'] = CalcMolFormula(mol)
    elif smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            metadata_dict['SMILES'] = Chem.MolToSmiles(mol)
            metadata_dict['INCHI'] = Chem.MolToInchi(mol)
            metadata_dict['INCHIKEY'] = Chem.MolToInchiKey(mol)
            metadata_dict['FORMULA'] = CalcMolFormula(mol)

    return metadata_dict

def mass_calculation(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information, including the SMILES string.
    :return: The updated metadata dictionary.

    This method performs mass calculations using the CDK library for the molecule specified by the SMILES string in the metadata dictionary. The method adds two new entries to the metadata
    * dictionary: 'EXACTMASS' and 'AVERAGEMASS', representing the exact mass and average mass of the molecule, respectively.

    Example usage:
    metadata = {
        "SMILES": "C1=CC=CC=C1"
    }
    mass_calculation(metadata)
    print(metadata)  # {'SMILES': 'C1=CC=CC=C1', 'EXACTMASS': '90.0399', 'AVERAGEMASS': '91.0938'}
    """
    SMILES = metadata_dict["SMILES"]

    if SMILES:
        # CDK exact mass
        try:
            mol = MolFromSmiles(SMILES)
            if mol is not None:
                metadata_dict['EXACTMASS'] = str(getMolExactMass(mol))
        except:
            return metadata_dict
        # CDK average mass
        try:
            mol = MolFromSmiles(SMILES)
            if mol is not None:
                metadata_dict['AVERAGEMASS'] = str(getMolNaturalMass(mol))
        except:
            return metadata_dict

    return metadata_dict

def mols_calculation(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata about a compound or molecule.

    :return: The updated metadata dictionary with additional derived properties.

    This method performs calculations on the given metadata to derive additional information about a compound or molecule. It first calls the 'mols_derivator' function to compute molecule
    * derivatives and updates the metadata dictionary accordingly. Then, it calls the 'mass_calculation' function to calculate the mass of the molecule and updates the metadata dictionary
    * again.

    The final updated metadata_dict is returned.
    """
    metadata_dict = mols_derivator(metadata_dict)
    metadata_dict = mass_calculation(metadata_dict)

    return metadata_dict