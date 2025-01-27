from GUI.utils.global_vars import parameters_dict
import deletion_report
import global_report


def calculate_unique_inchikeys(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC,
                               NEG_GC_insilico):
    """
    Calculate the number of unique INCHIKEYs in each DataFrame and update global_report.report_dict.

    Each key in the report_dict corresponds to a specific subset of the data.

    :param POS_LC: DataFrame containing positive LC experimental data.
    :param POS_LC_insilico: DataFrame containing positive LC in-silico data.
    :param POS_GC: DataFrame containing positive GC experimental data.
    :param POS_GC_insilico: DataFrame containing positive GC in-silico data.
    :param NEG_LC: DataFrame containing negative LC experimental data.
    :param NEG_LC_insilico: DataFrame containing negative LC in-silico data.
    :param NEG_GC: DataFrame containing negative GC experimental data.
    :param NEG_GC_insilico: DataFrame containing negative GC in-silico data.
    """

    # Fonction utilitaire pour compter les `INCHIKEY` uniques
    def count_unique_inchikeys(df):
        if 'INCHIKEY' in df.columns:  # Vérification si la colonne INCHIKEY existe
            return df['INCHIKEY'].nunique()
        return 0  # Retourne 0 si la colonne est absente

    # Mise à jour du dictionnaire report_dict avec les valeurs calculées
    global_report.report_dict["pos_lc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_LC)
    global_report.report_dict["neg_lc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_LC)
    global_report.report_dict["pos_lc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_LC_insilico)
    global_report.report_dict["neg_lc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_LC_insilico)
    global_report.report_dict["pos_gc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_GC)
    global_report.report_dict["neg_gc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_GC)
    global_report.report_dict["pos_gc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_GC_insilico)
    global_report.report_dict["neg_gc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_GC_insilico)


def calculate_spectrum_number(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC,
                              NEG_GC_insilico):
    """
    Calculate the number of spectra (rows) in each DataFrame and update global_report.report_dict.

    :param POS_LC: DataFrame containing positive LC experimental data.
    :param POS_LC_insilico: DataFrame containing positive LC in-silico data.
    :param POS_GC: DataFrame containing positive GC experimental data.
    :param POS_GC_insilico: DataFrame containing positive GC in-silico data.
    :param NEG_LC: DataFrame containing negative LC experimental data.
    :param NEG_LC_insilico: DataFrame containing negative LC in-silico data.
    :param NEG_GC: DataFrame containing negative GC experimental data.
    :param NEG_GC_insilico: DataFrame containing negative GC in-silico data.
    """

    # Mise à jour des clés dans le dictionnaire global `report_dict`
    global_report.report_dict["pos_lc_exp_spectrum_number"] = len(POS_LC)
    global_report.report_dict["neg_lc_exp_spectrum_number"] = len(NEG_LC)
    global_report.report_dict["pos_lc_insilico_spectrum_number"] = len(POS_LC_insilico)
    global_report.report_dict["neg_lc_insilico_spectrum_number"] = len(NEG_LC_insilico)
    global_report.report_dict["pos_gc_exp_spectrum_number"] = len(POS_GC)
    global_report.report_dict["neg_gc_exp_spectrum_number"] = len(NEG_GC)
    global_report.report_dict["pos_gc_insilico_spectrum_number"] = len(POS_GC_insilico)
    global_report.report_dict["neg_gc_insilico_spectrum_number"] = len(NEG_GC_insilico)


def format_parameters():
    parameters_string = f""" 
    ======================= FILES =======================
    INPUT_FILES:
        MSP:
    {"".join([f"\t\t\t{file}\n" for file in parameters_dict["input_directory"] if file.endswith(".msp")]) or "\t\t\t-- no file --\n"}
        JSON:
    {"".join([f"\t\t\t{file}\n" for file in parameters_dict["input_directory"] if file.endswith(".json")]) or "\t\t\t-- no file --\n"}
        CSV:
    {"".join([f"\t\t\t{file}\n" for file in parameters_dict["input_directory"] if file.endswith(".csv")]) or "\t\t\t-- no file --\n"}
        MGF:
    {"".join([f"\t\t\t{file}\n" for file in parameters_dict["input_directory"] if file.endswith(".mgf")]) or "\t\t\t-- no file --\n"}

    OUTPUT_DIRECTORY:
        {parameters_dict["output_directory"]}

    ===================== PARAMETERS =====================
    normalize_intensity: {"ON" if parameters_dict["normalize_intensity"] else "OFF"}
    remove_peak_above_precursormz: {"ON" if parameters_dict['remove_peak_above_precursormz'] else "OFF"}
    check_minimum_peak_requiered: {"ON" if parameters_dict['check_minimum_peak_requiered'] else "OFF"}
    reduce_peak_list: {"ON" if parameters_dict['reduce_peak_list'] else "OFF"}
    remove_spectrum_under_entropy_score: {"ON" if parameters_dict['remove_spectrum_under_entropy_score'] else "OFF"}
    keep_mz_in_range: {"ON" if parameters_dict['keep_mz_in_range'] else "OFF"}
    check_minimum_of_high_peaks_requiered: {"ON" if parameters_dict['check_minimum_of_high_peaks_requiered'] else "OFF"}
    
    reset_updates: {"YES" if parameters_dict['reset_updates'] else "NO"}
    
    """

    return parameters_string


def format_fitered_out():
    pass

def format_spectrum_numbers():
    pass

def format_unique_inchikeys():
    pass

def format_report():
    pass


def report(output_directory, POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico):
    calculate_spectrum_number(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico)
    calculate_unique_inchikeys(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico)
    format_report()
