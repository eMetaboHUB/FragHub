from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.Descriptors import ExactMolWt, MolWt
from rdkit import RDLogger, Chem

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
            mol = Chem.MolFromInchi(transforms['INCHI']) if 'InChI=' in inchi_smiles else Chem.MolFromSmiles(transforms['SMILES'])
            if mol is not None:
                try:
                    transforms['EXACTMASS'] = ExactMolWt(mol)
                    transforms['AVERAGEMASS'] = MolWt(mol)
                except:
                    return transforms

    return transforms


inchi_smiles = "COC(=O)C1=COC(OC2OC(CO)C(O)C(O)C2O)C2C1C(O)C=C2COC(=O)\C=C\c1ccc(O)cc1"
print(apply_transformations(inchi_smiles))
