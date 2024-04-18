import pandas as pd

def structure_traces_one(spectrum_list_TRACES):
    """
    Create a list of dictionaries containing structured traces from a given list of spectra.

    :param spectrum_list_TRACES: A list of spectra traces.
    :type spectrum_list_TRACES: list
    :return: A list of dictionaries containing structured traces.
    :rtype: list
    """
    TRACES_DICT_LIST = []

    for spectrums_traces in spectrum_list_TRACES:
        if spectrums_traces[0] is not None:
            TRACES_DICT_LIST.append({'old_spectrum': spectrums_traces[5],
                                     'no_smiles_no_inchi': spectrums_traces[1],
                                     'no_or_bad_precursormz_and_no_or_bad_addcut': spectrums_traces[2],
                                     'no_peaks_list': spectrums_traces[3],
                                     'minimum_peaks_not_requiered': spectrums_traces[4],
                                     })

    return TRACES_DICT_LIST

def structure_traces_two(spectrum_list_TRACES, DELETED_CONCATENATE_DF, profile_name):
    """
    :param spectrum_list_TRACES: A list containing dictionaries of spectrums. Each dictionary represents a spectrum and should have keys 'old_spectrum', 'no_smiles_no_inchi', 'no_or_bad_precursormz_and_no_or_bad_addcut', 'no_peaks_list', and 'minimum_peaks_not_requiered'. The value of 'old_spectrum' key should be a dictionary containing spectrum information.
    :param DELETED_CONCATENATE_DF: A DataFrame containing deleted spectrums.
    :param profile_name: A string representing the profile name.
    :return: None

    This method takes in a list of spectrums, merges them into a DataFrame, concatenates it with DELETED_CONCATENATE_DF, and saves the resulting DataFrame as an Excel file in the specified directory.
    """
    spectrum_list_TRACES_DF = pd.DataFrame()

    for spectrums in spectrum_list_TRACES:
        temp_spectrums = spectrums["old_spectrum"]
        temp_spectrums.update({"no_smiles_no_inchi": spectrums["no_smiles_no_inchi"],
                               "no_or_bad_precursormz_and_no_or_bad_addcut": spectrums["no_or_bad_precursormz_and_no_or_bad_addcut"],
                               "no_peaks_list": spectrums["no_peaks_list"],
                               "minimum_peaks_not_requiered": spectrums["minimum_peaks_not_requiered"],
                               "will_be_deleted_because_of_RDkit": False
                               })

        temp_df = pd.DataFrame(temp_spectrums, index=[0])

        spectrum_list_TRACES_DF = pd.concat([spectrum_list_TRACES_DF, temp_df])

    DF = pd.concat([spectrum_list_TRACES_DF, DELETED_CONCATENATE_DF])

    DF.to_excel(rf"../OUTPUT/{profile_name}/TRACES.xlsx", index=False)


