import matchms.filtering as msfilters
from standardizers import *
import concurrent.futures
import matchms.Fragments
import matchms.Metadata
import matchms.hashing
from tqdm import tqdm
from filters import *


def matchms_spectrum_to_str_msp(spectrum,file_name):
    """
    :param spectrum: The input spectrum to be converted to MSP format.
    :param file_name: The name of the file for the spectrum.
    :return: The spectrum converted to MSP format as a string.

    This method takes a spectrum object and a file name, and converts the spectrum object to MSP format. The MSP format is a tab-separated text format commonly used for representing mass
    * spectrometry data. The method first checks if the spectrum object is not None. If it is not None, it checks if the "filename" key is present in the metadata of the spectrum object
    *. If it is not present, it adds the file name to the output MSP string. Then, it iterates through the metadata items of the spectrum object and adds them to the output MSP string, excluding
    * keys that start with "num peaks", "num_peaks", or "peak_comments". After adding the metadata, it adds the number of peaks in the spectrum to the output MSP string. Finally, it iter
    *ates through the mz and intensity values of the peaks in the spectrum and adds them to the output MSP string in tab-separated format.

    Example usage:
    spectrum = ...
    file_name = "example.msp"
    msp_string = matchms_spectrum_to_str_msp(spectrum, file_name)
    print(msp_string)

    Output:
    FILENAME: example.msp
    METADATA_KEY1: value1
    METADATA_KEY2: value2
    ...
    NUM PEAKS: 100
    mz1         intensity1
    mz2         intensity2
    ...
    """
    if spectrum is not None:
        if "filename" not in  spectrum.metadata.keys():
            SPECTRUM = "FILENAME: " + file_name + "\n"
        else:
            SPECTRUM = ""
        for key, value in spectrum.metadata.items():
            if not (key.lower().startswith("num peaks") or key.lower().startswith("num_peaks") or key.lower().startswith("peak_comments")):
                KEY = re.sub(r'^[\W_]+|[\W_]+$', '', key.upper())
                SPECTRUM += f"{KEY}: {value}\n"

        SPECTRUM += f"NUM PEAKS: {len(spectrum.peaks)}\n"

        for mz, intensity in zip(spectrum.peaks.mz, spectrum.peaks.intensities):
            SPECTRUM += f"{mz}\t{intensity}\n".expandtabs(12)

        return SPECTRUM



def multithreaded_matchms(spectrum,file_name):
    """

    This method `multithreaded_matchms` takes two parameters: `spectrum` and `file_name`.

    :param spectrum: The spectrum object on which various metadata filters and peak filters will be applied.
    :param file_name: The name of the file associated with the spectrum.

    :return: The processed spectrum object with applied metadata and peak filters.

    The method proceeds with the following steps:

    1. Applies default metadata filters using `msfilters.default_filters`.

    2. Adds parent mass information to the spectrum using `msfilters.add_parent_mass`, estimating from the adduct if available.

    3. Converts the existing precursor mz field to a float using `msfilters.add_precursor_mz`.

    4. Adds retention time information to the spectrum using `msfilters.add_retention.add_retention_time`.

    5. Cleans the adduct information in the spectrum using `matchms.metadata_utils.clean_adduct`.

    6. Repairs the inchi, inchikey, and smiles fields in the spectrum using `msfilters.repair_inchi_inchikey_smiles`. For example, if the inchi field is present, it will be used to update
    * the inchikey field.

    7. Normalizes the intensities of the spectrum using `msfilters.normalize_intensities`.

    8. Selects peaks based on their relative intensity using `msfilters.select_by_relative_intensity`, with a threshold of 0.01 and a maximum intensity of 1.

    9. Selects peaks within a specified m/z range using `msfilters.select_by_mz`, with a minimum m/z value of 50 and a maximum m/z value of 2000.0.

    10. Reduces the number of peaks in the spectrum to a maximum of 500 using `msfilters.reduce_to_number_of_peaks`.

    11. Requires a minimum number of peaks in the spectrum using `msfilters.require_minimum_number_of_peaks`, with a minimum number of 3.

    12. Requires a minimum number of high peaks in the spectrum using `msfilters.require_minimum_of_high_peaks`, with a minimum number of peaks of 2 and an intensity percentage of 5.0.

    13. Converts the spectrum object to a string in MSP format using `matchms_spectrum_to_str_msp`.

    14. Harmonizes the field names in the spectrum using `harmonize_fields_names`.

    15. Harmonizes the field values in the spectrum using `harmonize_fields_values`.

    16. Returns the processed and filtered spectrum object.

    """
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

    spectrum = matchms_spectrum_to_str_msp(spectrum,file_name)

    spectrum = harmonize_fields_names(spectrum)

    spectrum = harmonize_fields_values(spectrum)

    spectrum = remove_peaks_above_precursor_mz(spectrum)

    return spectrum

def matchms_processing(spectrum_list,file_name):
    """
    Process the given spectrum list using multithreading.

    :param spectrum_list: The list of spectra to process.
    :param file_name: The name of the file being processed.
    :return: The list of different worker executions.
    """
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(multithreaded_matchms, spectrum_list, [file_name for i in range(len(spectrum_list))]), total=len(spectrum_list), unit=" spectrums", colour="green", desc="\t  processing"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


