import pandas as pd

def structure_traces_one(spectrum_list_TRACES, profile_name):
    """
    Create a list of dictionaries containing structured traces from a given list of spectra.

    :param spectrum_list_TRACES: A list of spectra traces.
    :type spectrum_list_TRACES: list
    :return: A list of dictionaries containing structured traces.
    :rtype: list
    """
    spectrum_list_TRACES_DF = pd.DataFrame()

    for spectrums in spectrum_list_TRACES:
        if spectrums[0] is None:
            temp_spectrums = spectrums[5]
            temp_spectrums.update({"no_smiles_no_inchi": spectrums[1],
                                   "no_or_bad_precursormz_and_no_or_bad_addcut": spectrums[2],
                                   "no_peaks_list": spectrums[3],
                                   "peaks_filters_not_valid": spectrums[4],
                                   "will_be_deleted_because_of_RDkit": False
                                   })

            # Used list-comprehension to ensure every value is a list
            temp_spectrums = {k: [v] for k, v in temp_spectrums.items()}

            temp_df = pd.DataFrame(temp_spectrums)  # Removed index, as we can use default index

            spectrum_list_TRACES_DF = pd.concat([spectrum_list_TRACES_DF, temp_df])

    spectrum_list_TRACES_DF.to_excel(rf"../OUTPUT/{profile_name}/TRACKER_1.xlsx", index=False)

    return spectrum_list_TRACES_DF

def structure_traces_two(spectrum_list_TRACES_DF, DELETED_CONCATENATE_DF, profile_name):
    """
    :param spectrum_list_TRACES: A list containing dictionaries of spectrums. Each dictionary represents a spectrum and should have keys 'old_spectrum', 'no_smiles_no_inchi', 'no_or_bad_precursormz_and_no_or_bad_addcut', 'no_peaks_list', and 'minimum_peaks_not_requiered'. The value of 'old_spectrum' key should be a dictionary containing spectrum information.
    :param DELETED_CONCATENATE_DF: A DataFrame containing deleted spectrums.
    :param profile_name: A string representing the profile name.
    :return: None

    This method takes in a list of spectrums, merges them into a DataFrame, concatenates it with DELETED_CONCATENATE_DF, and saves the resulting DataFrame as an Excel file in the specified directory.
    """
    DF = pd.concat([spectrum_list_TRACES_DF, DELETED_CONCATENATE_DF])

    DF.to_excel(rf"../OUTPUT/{profile_name}/TRACKER_2.xlsx", index=False)


