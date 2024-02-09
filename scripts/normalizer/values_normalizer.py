from .normalize_instruments_and_resolution import *
from .missing_precursormz_re_calculation import *
from .delete_no_smiles_no_inchi import *
from .normalize_retentiontime import *
from .repair_mol_descriptors import *
from .normalize_ionization import *
from .normalize_predicted import *
from .normalize_ms_level import *
from .normalize_empties import *
from .normalize_ionmode import *
from .normalize_adduct import *

def normalize_values(metadata_dict):
    """
    :param metadata_dict: A dictionary containing metadata information.
    :return: The normalized metadata dictionary.
    """
    metadata_dict = normalize_empties(metadata_dict)

    metadata_dict = repair_mol_descriptors(metadata_dict)

    metadata_dict = delete_no_smiles_no_inchi(metadata_dict)

    if metadata_dict:
        metadata_dict = normalize_adduct(metadata_dict)
        metadata_dict = missing_precursormz_re_calculation(metadata_dict)
        metadata_dict = normalize_ionmode(metadata_dict)
        metadata_dict = normalize_ms_level(metadata_dict)
        metadata_dict = normalize_predicted(metadata_dict)
        metadata_dict = normalize_retentiontime(metadata_dict)
        metadata_dict = normalize_ionization(metadata_dict)
        metadata_dict = normalize_instruments_and_resolution(metadata_dict)

    return metadata_dict
