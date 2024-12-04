import re

global smiles_pattern
smiles_pattern = re.compile(r"\b([^J][0-9BCOHNSOPrIFla@+\-\[\]\(\)\\\/%=#$]{6,})\b", flags=re.IGNORECASE) # Match smiles

global inchi_pattern
inchi_pattern = re.compile(r"(InChI=.*|\/[0-9A-Z]*\/)", flags=re.IGNORECASE) # Match inchi

global inchikey_pattern
inchikey_pattern = re.compile(r"([A-Z]{14}-[A-Z]{10}-[NO])", flags=re.IGNORECASE) # Match inchikey

global repair_inchi_pattern
repair_inchi_pattern = re.compile(r"(inchi=)", flags=re.IGNORECASE)

def repair_inchi(metadata_dict):
    """
    This function 'repair_inchi' is used to fix the 'INCHI' value in the given metadata dictionary if it exists.

    The INCHI string is fixed by replacing the leading pattern determined by 'repair_inchi_pattern' regex with "InChI=".

    :param metadata_dict: A dictionary containing metadata of a molecule. The 'INCHI' key in this dictionary is updated.

    :return: The same metadata dictionary with the updated 'INCHI' key if it was present initially.

    Note: If the 'INCHI' key is not in the dictionary or its value was None, then the dictionary is returned as it is, without any modification.
    """

    # Extraction of InChI string from metadata dictionary
    inchi = metadata_dict['INCHI']

    # If InChI string is present, the function modifies it.
    if inchi:
        # Replacing the pattern with "InChI=" in the existing InChI string
        inchi = re.sub(repair_inchi_pattern, "InChI=", inchi)

        # Updating the InChI key in metadata dictionary with the modified InChI string
        metadata_dict['INCHI'] = inchi

    # Metadata dictionary is returned (modified or unmodified)
    return metadata_dict

def repair_mol_descriptors(metadata_dict):
    """
    Repairs molecular descriptors in a dictionary containing SMILES, InChI, and InChIKey.
    Check if molecular descriptors are in the dedicated field.
    """
    smiles = metadata_dict['SMILES']
    inchi = metadata_dict['INCHI']
    inchikey = metadata_dict['INCHIKEY']
    name = metadata_dict['NAME']
    comment = metadata_dict['COMMENT']

    # ============= CASES =============
    # INCHIKEY
    inchikey_in_inchikey = re.search(inchikey_pattern, inchikey)
    inchikey_in_smiles = re.search(inchikey_pattern, smiles)
    inchikey_in_inchi = re.search(inchikey_pattern, inchi)
    inchikey_in_name = re.search(inchikey_pattern, name)
    inchikey_in_comment = re.search(inchikey_pattern, comment)

    if inchikey_in_inchikey:
        metadata_dict['INCHIKEY'] = inchikey_in_inchikey.group(1)
    if inchikey_in_smiles:
        metadata_dict['INCHIKEY'] = inchikey_in_smiles.group(1)
        metadata_dict['SMILES'] = smiles.replace(inchikey_in_smiles.group(1), '')
    if inchikey_in_inchi:
        metadata_dict['INCHIKEY'] = inchikey_in_inchi.group(1)
        metadata_dict['INCHI'] = inchi.replace(inchikey_in_inchi.group(1), '')
    if inchikey_in_name:
        metadata_dict['INCHIKEY'] = inchikey_in_name.group(1)
        metadata_dict['NAME'] = name.replace(inchikey_in_name.group(1), '')
    if inchikey_in_comment:
        metadata_dict['INCHIKEY'] = inchikey_in_comment.group(1)
        metadata_dict['COMMENT'] = comment.replace(inchikey_in_comment.group(1), '')

    # INCHI
    inchi_in_inchi = re.search(inchi_pattern, inchi)
    inchi_in_smiles = re.search(inchi_pattern, smiles)
    inchi_in_inchikey = re.search(inchi_pattern, inchikey)
    inchi_in_name = re.search(inchi_pattern, name)
    inchi_in_comment = re.search(inchi_pattern, comment)

    if inchi_in_inchi:
        metadata_dict['INCHI'] = inchi_in_inchi.group(1)
    if inchi_in_smiles:
        metadata_dict['INCHI'] = inchi_in_smiles.group(1)
        metadata_dict['SMILES'] = smiles.replace(inchi_in_smiles.group(1), '')
    if inchi_in_inchikey:
        metadata_dict['INCHI'] = inchi_in_inchikey.group(1)
        metadata_dict['INCHIKEY'] = inchikey.replace(inchi_in_inchikey.group(1), '')
    if inchi_in_name:
        metadata_dict['INCHI'] = inchi_in_name.group(1)
        metadata_dict['NAME'] = name.replace(inchi_in_name.group(1), '')
    if inchi_in_comment:
        metadata_dict['INCHI'] = inchi_in_comment.group(1)
        metadata_dict['COMMENT'] = comment.replace(inchi_in_comment.group(1), '')

    # SMILES
    smiles_in_smiles = re.search(smiles_pattern, smiles)
    smiles_in_inchi = re.search(smiles_pattern, inchi)
    smiles_in_inchikey = re.search(smiles_pattern, inchikey)
    smiles_in_name = re.search(smiles_pattern, name)
    smiles_in_comment = re.search(smiles_pattern, comment)

    if smiles_in_smiles:
        metadata_dict['SMILES'] = smiles_in_smiles.group(1)
    if smiles_in_inchi:
        metadata_dict['SMILES'] = smiles_in_inchi.group(1)
        metadata_dict['INCHI'] = inchi.replace(smiles_in_inchi.group(1), '')
    if smiles_in_inchikey:
        metadata_dict['SMILES'] = smiles_in_inchikey.group(1)
        metadata_dict['INCHIKEY'] = inchikey.replace(smiles_in_inchikey.group(1), '')
    if smiles_in_name:
        metadata_dict['SMILES'] = smiles_in_name.group(1)
        metadata_dict['NAME'] = name.replace(smiles_in_name.group(1), '')
    if smiles_in_comment:
        metadata_dict['SMILES'] = smiles_in_comment.group(1)
        metadata_dict['COMMENT'] = comment.replace(smiles_in_comment.group(1), '')

    # Again call the repair InChI function to repair and update the InChI
    metadata_dict = repair_inchi(metadata_dict)

    # Return the updated metadata dictionary
    return metadata_dict
