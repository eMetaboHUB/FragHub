from matchms.logging_functions import set_matchms_logger_level, add_logging_to_file
from matchms.importing import load_from_msp
from matchms.exporting import save_as_msp
from tqdm.notebook import tqdm as tqdm
import matchms.filtering as msfilters
import matchms.metadata_utils
import matchms.Fragments
import matchms.Metadata
import matchms.hashing
import pandas as pd
import re

import concurrent.futures

clean_msp_list = []

def harmonize_fields_names(spectrum):
    df = pd.read_csv("./data/colomnsNames.csv", sep=";")
    dict = df.to_dict()

    # re-structure dict
    for column, elements in dict.items():
        dict[column] = [element for key, element in elements.items() if (element != column) and (str(element) != 'nan')]

    # replace non_normalized fields by normalized fields
    for column, lists in dict.items():
        for non_normalize in lists:
            spectrum = re.sub(non_normalize, column, spectrum)

    return spectrum


def multithreaded_matchms(spectrum):
    # apply_metadata_filters
    spectrum = msfilters.default_filters(spectrum)
    spectrum = msfilters.derive_adduct_from_name(spectrum)
    spectrum = msfilters.add_parent_mass(spectrum, estimate_from_adduct=True)
    spectrum = msfilters.add_retention.add_retention_time(spectrum)
    spectrum = msfilters.set_ionmode_na_when_missing(spectrum)
    spectrum = matchms.metadata_utils.clean_adduct(spectrum)
    spectrum = msfilters.repair_inchi_inchikey_smiles(spectrum)
    spectrum = msfilters.derive_inchi_from_smiles(spectrum)
    spectrum = msfilters.derive_smiles_from_inchi(spectrum)
    spectrum = msfilters.derive_inchikey_from_inchi(spectrum)
    spectrum = msfilters.harmonize_undefined_smiles(spectrum)
    spectrum = msfilters.harmonize_undefined_inchi(spectrum)
    spectrum = msfilters.harmonize_undefined_inchikey(spectrum)
    spectrum = msfilters.add_precursor_mz(spectrum)

    # normalize_and_filter_peaks
    spectrum = msfilters.normalize_intensities(spectrum)
    spectrum = msfilters.select_by_relative_intensity(spectrum, 0.01, 1)
    spectrum = msfilters.select_by_mz(spectrum, mz_from=50, mz_to=2000.0)
    spectrum = msfilters.reduce_to_number_of_peaks(spectrum, n_max=500)
    spectrum = msfilters.require_minimum_number_of_peaks(spectrum, n_required=3)

    # Replace non-normalized fields names by normalized.
    spectrum = harmonize_fields_names(spectrum)

    return spectrum


def app_msp(spectrum):
    clean_msp_list.append(spectrum)

def matchms_treatment(spectrum_list):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        future_proc = {executor.submit(multithreaded_matchms, spectrum): spectrum for spectrum in spectrum_list}
        for future in concurrent.futures.as_completed(future_proc):
            app_msp(future.result())
