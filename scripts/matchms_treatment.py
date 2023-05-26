from matchms.logging_functions import set_matchms_logger_level, add_logging_to_file
from matchms.importing import load_from_msp
from matchms.exporting import save_as_msp
from tqdm.notebook import tqdm as tqdm
import matchms.filtering as msfilters
import matchms.metadata_utils
from matchms import Spectrum
import concurrent.futures
import matchms.Fragments
import matchms.Metadata
import matchms.hashing
import threading
import re

file_lock = threading.Lock()


def multithreaded_matchms(spectrum):
    thread_num = threading.get_ident()

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

    # Write to file with lock synchronization.
    # In this way, each thread will have to wait until the lock is available before it can enter the critical section (writing to the file).
    # This ensures that only one thread can write to the file at a time, avoiding conflicts.
    with file_lock:
        save_as_msp(spectrum, "./temp/"+str(thread_num)+"_temp.msp")

def matchms_treatment(spectrum_list):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        results = executor.map(multithreaded_matchms, spectrum_list)

