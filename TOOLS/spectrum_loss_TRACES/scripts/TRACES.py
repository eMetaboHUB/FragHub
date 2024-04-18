
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