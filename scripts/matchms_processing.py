import matchms.filtering as msfilters
import matchms.metadata_utils
from standardizers import *
import concurrent.futures
import matchms.Fragments
import matchms.Metadata
import matchms.hashing
from tqdm import tqdm
from filters import *


def matchms_spectrum_to_str_msp(spectrum,file_name):
    """
    :param spectrum: The spectrum to be converted to a string in the MSP format.
    :param file_name: The name of the file from which the spectrum originated.
    :return: The spectrum in string format following the MSP format.

    This method takes a spectrum object and converts it into a string formatted according to the MSP (Mass Spectral Peak) format. The spectrum object should have the following attributes
    *:

    - `metadata`: A dictionary containing metadata information for the spectrum.
    - `peaks`: A list-like object containing the mass-to-charge ratio (m/z) values and intensities of the peaks in the spectrum.

    The method checks if the spectrum object is not None, and if the key "filename" is not present in the metadata. If the key is not present, it adds a line to the string indicating the
    * file name.

    Next, the method iterates through the metadata items and excludes any items with keys that start with "num peaks", "num_peaks", or "peak_comments". It appends the remaining metadata
    * items to the string.

    After that, the method appends the number of peaks in the spectrum to the string.

    Finally, the method appends each m/z value and intensity pair from the spectrum peaks to the string, with the values formatted and aligned for readability.

    The resulting string in the MSP format is returned.
    """
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
    """
    :param spectrum: The spectrum object to be processed.
    :param file_name: The name of the file containing the spectrum.
    :return: The processed spectrum object.

    This method applies a series of filters and transformations to the input spectrum object to preprocess it. The processed spectrum is then returned.

    The processing steps include:
    - Applying metadata filters using the msfilters.default_filters() function.
    - Adding parent mass information using msfilters.add_parent_mass() function.
    - Adding precursor m/z information using msfilters.add_precursor_mz() function.
    - Adding retention time information using msfilters.add_retention_time() function.
    - Cleaning adduct information using matchms.metadata_utils.clean_adduct() function.
    - Repairing InChi, InChiKey, and SMILES information using msfilters.repair_inchi_inchikey_smiles() function.
    - Normalizing intensities using msfilters.normalize_intensities() function.
    - Selecting peaks based on relative intensity using msfilters.select_by_relative_intensity() function.
    - Selecting peaks based on m/z range using msfilters.select_by_mz() function.
    - Reducing the number of peaks to a maximum number using msfilters.reduce_to_number_of_peaks() function.
    - Requiring a minimum number of peaks using msfilters.require_minimum_number_of_peaks() function.
    - Requiring a minimum number of high peaks based on intensity using msfilters.require_minimum_of_high_peaks() function.
    - Converting the spectrum object to string MSP format using the matchms_spectrum_to_str_msp() function.
    - Harmonizing the field names of the spectrum object using the harmonize_fields_names() function.
    - Harmonizing the field values of the spectrum object using the harmonize_fields_values() function.

    The processed spectrum object is then returned.
    """
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
    """
    Perform multithreaded processing of spectrum_list using matchms.

    :param spectrum_list: List of spectra to process.
    :type spectrum_list: list
    :param file_name: Name of the file being processed.
    :type file_name: str
    :return: List of different worker executions.
    :rtype: list
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(multithreaded_matchms, spectrum_list, [file_name for i in range(len(spectrum_list))]), total=len(spectrum_list), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


