import re

global smiles_pattern
smiles_pattern = re.compile(r"[^J][a-z0-9@+\-\[\]\(\)\\\/%=#$]{6,}", flags=re.IGNORECASE) # Match smiles

global inchi_pattern
inchi_pattern = re.compile(r"InChI=.*|\/[0-9A-Z]*\/", flags=re.IGNORECASE) # Match inchi

global inchikey_pattern
inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])|([A-Z]{14})", flags=re.IGNORECASE) # Match inchikey or short inchikey

global repair_inchi_pattern
repair_inchi_pattern = re.compile(r"^(inchi=)?", flags=re.IGNORECASE)

def repair_inchi(metadata_dict):
    """
    :param metadata_dict: dictionary containing metadata
    :return: modified metadata dictionary with 'INCHI' key updated

    """
    inchi = metadata_dict['INCHI']

    if inchi:
        inchi = re.sub(repair_inchi_pattern,"InChI=",inchi)
        metadata_dict['INCHI'] = inchi

        return metadata_dict

    return metadata_dict

def repair_mol_descriptors(metadata_dict):
    """
        repair_mol_descriptors(metadata_dict)

        Repairs molecular descriptors in a dictionary containing SMILES, InChI, and InChIKey.
        Check if molecular descriptors are in the dedicated field.

        :param metadata_dict: A dictionary containing molecular descriptors.
        :return: The repaired dictionary with updated molecular descriptors.

        Example Usage:

        metadata_dict = {
            'SMILES': 'CC(=O)O',
            'INCHI': 'InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)'
            'INCHIKEY': 'QTBSBXVTEAMEQO-UHFFFAOYSA-N'
        }

        repaired_dict = repair_mol_descriptors(metadata_dict)

        The repaired_dict will contain the repaired molecular descriptors.
    """
    smiles = metadata_dict['SMILES']
    inchi = metadata_dict['INCHI']
    inchikey = metadata_dict['INCHIKEY']

    if re.search(smiles_pattern, smiles) and re.search(inchi_pattern, inchi) and re.search(inchikey_pattern, inchikey):
        metadata_dict = repair_inchi(metadata_dict)
        return metadata_dict

    # SMILES
    if re.search(smiles_pattern, inchi):
        if not re.search(inchi_pattern, inchi) and not re.search(inchikey_pattern, inchi):
            metadata_dict['SMILES'] = inchi
            metadata_dict['INCHI'] = ''
    if re.search(smiles_pattern, inchikey):
        if not re.search(inchi_pattern, inchikey) and not re.search(inchikey_pattern, inchikey):
            metadata_dict['SMILES'] = inchikey
            metadata_dict['INCHIKEY'] = ''
    # INCHI
    if re.search(inchi_pattern, smiles):
        metadata_dict['INCHI'] = smiles
        metadata_dict['INCHI'] = ''
    if re.search(inchi_pattern, inchikey):
        metadata_dict['INCHI'] = inchikey
        metadata_dict['INCHIKEY'] = ''
    # INCHIKEY
    if re.search(inchikey_pattern, smiles):
        metadata_dict['INCHIKEY'] = smiles
        metadata_dict['SMILES'] = ''
    if re.search(inchikey_pattern, inchi):
        metadata_dict['INCHIKEY'] = inchi
        metadata_dict['INCHI'] = ''

    metadata_dict = repair_inchi(metadata_dict)

    return metadata_dict