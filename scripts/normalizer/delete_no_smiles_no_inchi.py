import deletion_report

def delete_no_smiles_no_inchi_no_inchikey(metadata_dict):
    """
    This function deletes entries from the provided `metadata_dict` dictionary
    if both 'SMILES' and 'INCHI' keys are not present in it (i.e., have NaN values).

    :param metadata_dict: A dictionary containing metadata about chemical compounds. It usually
                          contains keys such as 'SMILES', 'INCHI', etc. Sometimes these keys
                          may have NaN values indicating their absence.
    :type metadata_dict: dict

    :return: If both 'SMILES' and 'INCHI' keys in the input dictionary have NaN values, None is returned.
             Otherwise, the original dictionary is returned, i.e., no entries are deleted.
    :rtype: dict or None
    """
    # Check if both 'SMILES' and 'INCHI' keys in the dictionary do not exist (have NaN values).
    if not metadata_dict["SMILES"] and not metadata_dict["INCHI"] and not metadata_dict["INCHIKEY"]:
        # If both keys do not exist, return None. This effectively deletes the entries from a
        # higher-level context as the returned None may not be added back to a collection of metadata dictionaries.
        metadata_dict['DELETION_REASON'] = "spectrum deleted because because it has neither inchi nor smiles nor inchikey"
        deletion_report.deleted_spectrum_list.append(metadata_dict)
        deletion_report.no_smiles_no_inchi_no_inchikey += 1
        return None
    else:
        # If either 'SMILES' or 'INCHI' key exists (the values are not NaN), return the original dictionary.
        return metadata_dict