from matchms.logging_functions import set_matchms_logger_level
from matchms.importing import load_from_msp
import matchms.filtering as msfilters
from matchms.exporting import *
import matchms.metadata_utils
from msp_utilities import *
import concurrent.futures
import matchms.Fragments
import matchms.Metadata
import matchms.hashing
from tqdm import tqdm
import os
import re

def matchms_spectrum_to_str_msp(spectrum,file_name):
    if spectrum is not None:
        if "filename" not in  spectrum.metadata.keys():
            SPECTRUM = "FILENAME: " + file_name + "\n"
        else:
            SPECTRUM = ""
        for key, value in spectrum.metadata.items():
            if not (key.lower().startswith("num peaks") or key.lower().startswith("num_peaks") or key.lower().startswith("peak_comments")):
                SPECTRUM += f"{key.upper()}: {value}\n"

        SPECTRUM += f"NUM PEAKS: {len(spectrum.peaks)}\n"

        for mz, intensity in zip(spectrum.peaks.mz, spectrum.peaks.intensities):
            SPECTRUM += f"{mz}\t{intensity}\n".expandtabs(12)

        return SPECTRUM



def multithreaded_matchms(spectrum,file_name):

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

    # normalize_and_filter_peaks
    spectrum = msfilters.normalize_intensities(spectrum)
    spectrum = msfilters.select_by_relative_intensity(spectrum, 0.01, 1)
    spectrum = msfilters.select_by_mz(spectrum, mz_from=50, mz_to=2000.0)
    spectrum = msfilters.reduce_to_number_of_peaks(spectrum, n_max=500)
    spectrum = msfilters.require_minimum_number_of_peaks(spectrum, n_required=3)

    spectrum = matchms_spectrum_to_str_msp(spectrum,file_name)

    spectrum = harmonize_fields_names(spectrum)
    
    spectrum = harmonize_fields_values(spectrum)

    return spectrum

def matchms_processing(spectrum_list,file_name):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(multithreaded_matchms, spectrum_list, [file_name for i in range(len(spectrum_list))]), total=len(spectrum_list), unit="spectrums", colour="green"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


