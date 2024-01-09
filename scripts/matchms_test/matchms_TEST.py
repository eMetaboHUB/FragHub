from matchms.importing import load_from_msp
from matchms.exporting import save_as_msp
import matchms.filtering as msfilters
import time
import os

start_time = time.time()

def matchms_filters(spectrum):
    # apply_metadata_filters
    spectrum = msfilters.default_filters(spectrum)
    spectrum = msfilters.add_parent_mass(spectrum, estimate_from_adduct=True)
    spectrum = msfilters.add_precursor_mz(spectrum) # ne fait que convertir le champ déjà existant en float, ne calcule rien
    spectrum = msfilters.add_retention_time(spectrum)
    spectrum = msfilters.clean_adduct(spectrum)
    spectrum = msfilters.repair_inchi_inchikey_smiles(spectrum) # example: si inchi dans champ inchikey

    # normalize_and_filter_peaks
    spectrum = msfilters.normalize_intensities(spectrum)
    spectrum = msfilters.select_by_relative_intensity(spectrum, 0.01, 1)
    spectrum = msfilters.select_by_mz(spectrum, mz_from=50, mz_to=2000.0)
    spectrum = msfilters.reduce_to_number_of_peaks(spectrum, n_max=500)
    spectrum = msfilters.require_minimum_number_of_peaks(spectrum, n_required=3)
    spectrum = msfilters.require_minimum_number_of_high_peaks(spectrum, no_peaks=2 ,intensity_percent=5.0)

    return spectrum

msp_dir = r'../../INPUT/MSP'

for files in os.listdir(msp_dir):
    if files.endswith('.msp'):
        msp_path = os.path.join(msp_dir,files)
        spectrum_list = list(load_from_msp(msp_path))

        # Apply filters to clean and enhance each spectrum
        spectrums_final = []
        for spectrum in spectrum_list:
            # Apply default filter to standardize ion mode, correct charge and more.
            # Default filter is fully explained at https://matchms.readthedocs.io/en/latest/api/matchms.filtering.html .
            spectrum = matchms_filters(spectrum)
            spectrums_final.append(spectrum)

        spectrums_final = [spectrum for spectrum in spectrums_final if spectrum is not None]

        save_as_msp(spectrums_final,r"C:\Users\Axel\PycharmProjects\msp_v3\scripts\matchms_test\ALL_GNPS.msp")

print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))