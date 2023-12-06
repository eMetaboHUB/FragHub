import matchms.filtering as msfilters
import matchms.metadata_utils
from standardizers import *
import concurrent.futures
import matchms.Fragments
import matchms.Metadata
import matchms.hashing
from tqdm import tqdm


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
    spectrum = msfilters.add_parent_mass(spectrum, estimate_from_adduct=True)
    spectrum = msfilters.add_precursor_mz(spectrum) # ne fait que convertir le champ déjà existant en float, ne calcule rien
    spectrum = msfilters.add_retention.add_retention_time(spectrum)
    spectrum = matchms.metadata_utils.clean_adduct(spectrum)
    spectrum = msfilters.repair_inchi_inchikey_smiles(spectrum) # example: si inchi dans champ inchikey

    # normalize_and_filter_peaks
    spectrum = msfilters.normalize_intensities(spectrum)
    spectrum = msfilters.select_by_relative_intensity(spectrum, 0.01, 1)
    spectrum = msfilters.select_by_mz(spectrum, mz_from=50, mz_to=2000.0)
    spectrum = msfilters.reduce_to_number_of_peaks(spectrum, n_max=500)
    spectrum = msfilters.require_minimum_number_of_peaks(spectrum, n_required=3)
    spectrum = msfilters.require_minimum_of_high_peaks(no_peaks=2 ,intensity_percent=5.0)

    spectrum = matchms_spectrum_to_str_msp(spectrum,file_name)

    spectrum = harmonize_fields_names(spectrum)

    spectrum = harmonize_fields_values(spectrum)

    return spectrum

def matchms_processing(spectrum_list,file_name):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(multithreaded_matchms, spectrum_list, [file_name for i in range(len(spectrum_list))]), total=len(spectrum_list), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


