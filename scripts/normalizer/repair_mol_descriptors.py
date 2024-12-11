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

    # Assigning the values of SMILES, InChI, and InChIKey into respective variables
    smiles = metadata_dict['SMILES']
    inchi = metadata_dict['INCHI']
    inchikey = metadata_dict['INCHIKEY']

    # If the respective patterns match the respective values, repair the InChI and return the updated dict
    if re.search(smiles_pattern, smiles) and re.search(inchi_pattern, inchi) and re.search(inchikey_pattern, inchikey):
        metadata_dict = repair_inchi(metadata_dict)
        return metadata_dict

    # If the patterns of InChI and InChIKey do not match but SMILES pattern matches InChI,
    # then update InChI value to SMILES and set the InChI to blank
    if re.search(smiles_pattern, inchi):
        if not re.search(inchi_pattern, inchi) and not re.search(inchikey_pattern, inchi):
            metadata_dict['SMILES'] = inchi
            metadata_dict['INCHI'] = ''

    # If the patterns of InChi and InChiKey do not match but SMILES pattern matches InChIKey,
    # update the InChIKey value to SMILES and set the InChIKey to blank
    if re.search(smiles_pattern, inchikey):
        if not re.search(inchi_pattern, inchikey) and not re.search(inchikey_pattern, inchikey):
            metadata_dict['SMILES'] = inchikey
            metadata_dict['INCHIKEY'] = ''

    # If SMILES matches InChI pattern, set the InChI to SMILES and set SMILES to blank
    if re.search(inchi_pattern, smiles):
        metadata_dict['INCHI'] = smiles
        metadata_dict['SMILES'] = ''

    # If InChIKey matches InChI pattern, set the InChI to InChIKey and set the InChiKey to blank
    if re.search(inchi_pattern, inchikey):
        metadata_dict['INCHI'] = inchikey
        metadata_dict['INCHIKEY'] = ''

    # If InChI matches InChiKey pattern, set the InChiKey to InChI and set InChI to blank
    if re.search(inchikey_pattern, inchi):
        metadata_dict['INCHIKEY'] = inchi
        metadata_dict['INCHI'] = ''

    if re.search(inchikey_pattern, smiles):
        metadata_dict['INCHIKEY'] = smiles
        metadata_dict['SMILES'] = ''

    # Again call the repair InChI function to repair and update the InChI
    metadata_dict = repair_inchi(metadata_dict)

    # Return the updated metadata dictionary
    return metadata_dict
