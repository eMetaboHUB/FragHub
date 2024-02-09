
def delete_no_smiles_no_inchi(metadata_dict):
    """
    Delete entries from the given metadata dictionary if both 'SMILES' and 'INCHI' keys have NaN values.

    :param metadata_dict: A dictionary containing metadata information.
    :type metadata_dict: dict
    :return: The updated metadata dictionary.
    :rtype: dict
    """
    if not metadata_dict["SMILES"] and not metadata_dict["INCHI"]:
        return None
    else:
        return metadata_dict