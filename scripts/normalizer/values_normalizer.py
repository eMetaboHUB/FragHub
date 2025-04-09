from .normalize_instruments_and_resolution import *
from .missing_precursormz_re_calculation import *
from .delete_no_smiles_no_inchi import *
from .normalize_retentiontime import *
from .repair_mol_descriptors import *
from .check_for_bad_adduct import *
from .normalize_ionization import *
from .normalize_predicted import *
from .normalize_ms_level import *
from .normalize_empties import *
from .normalize_ionmode import *
from .normalize_adduct import *

def normalize_values(metadata_dict):
    """
    This function takes in a metadata dictionary and applies numerous normalization functions on it to standardize its values.

    :param metadata_dict: A dictionary containing metadata information.
    :return: The normalized metadata dictionary.
    """
    # Replace any missing or 'NaN' values with appropriate values
    metadata_dict = normalize_empties(metadata_dict)

    # Repair molecular descriptors in the metadata
    metadata_dict = repair_mol_descriptors(metadata_dict)

    # If both 'SMILES' and 'INCHI' keys in the dictionary
    # do not exist (have NaN values), we delete the dictionary.
    metadata_dict = delete_no_smiles_no_inchi_no_inchikey(metadata_dict)

    # If after the above operations the metadata_dict is not empty continue with the normalization
    if metadata_dict:
        # Normalize Ionization in the metadata
        metadata_dict = normalize_ionization(metadata_dict)

        # Normalize the instruments and resolution data in the metadata
        metadata_dict = normalize_instruments_and_resolution(metadata_dict)

        # Normalize adduct data in the metadata
        metadata_dict = normalize_adduct(metadata_dict)

        # If 'PRECURSORMZ' key does not exist in the metadata_dict
        # or no recalculation is needed, perform recalculation
        metadata_dict = missing_precursormz_re_calculation(metadata_dict)

        # Normalize the ion mode in the metadata, from long form to short standardized form
        metadata_dict = normalize_ionmode(metadata_dict)

        # ckeck if adduct in pos is really pos (exemple)
        print(metadata_dict, "\n\n")
        metadata_dict = check_for_bad_adduct(metadata_dict)

        if metadata_dict:

            # Normalize MS level
            metadata_dict = normalize_ms_level(metadata_dict)

            # Normalize the predicted value in the metadata
            metadata_dict = normalize_predicted(metadata_dict)

            # Normalize Retention Time in the metadata which can be represented in different units
            metadata_dict = normalize_retentiontime(metadata_dict)


    return metadata_dict
