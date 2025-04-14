from scripts.GUI.utils.global_vars import parameters_dict
from datetime import datetime
import scripts.deletion_report
import scripts.global_report
import os

def calculate_unique_inchikeys(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df):
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
    scripts.global_report.report_dict["pos_lc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_LC_df)
    scripts.global_report.report_dict["neg_lc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_LC_df)
    scripts.global_report.report_dict["pos_lc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_LC_In_Silico_df)
    scripts.global_report.report_dict["neg_lc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_LC_In_Silico_df)
    scripts.global_report.report_dict["pos_gc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_GC_df)
    scripts.global_report.report_dict["neg_gc_exp_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_GC_df)
    scripts.global_report.report_dict["pos_gc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(POS_GC_In_Silico_df)
    scripts.global_report.report_dict["neg_gc_insilico_spectrum_unique_inchikey"] = count_unique_inchikeys(NEG_GC_In_Silico_df)


def calculate_spectrum_number(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df):
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
    scripts.global_report.report_dict["pos_lc_exp_spectrum_number"] = len(POS_LC_df)
    scripts.global_report.report_dict["neg_lc_exp_spectrum_number"] = len(NEG_LC_df)
    scripts.global_report.report_dict["pos_lc_insilico_spectrum_number"] = len(POS_LC_In_Silico_df)
    scripts.global_report.report_dict["neg_lc_insilico_spectrum_number"] = len(NEG_LC_In_Silico_df)
    scripts.global_report.report_dict["pos_gc_exp_spectrum_number"] = len(POS_GC_df)
    scripts.global_report.report_dict["neg_gc_exp_spectrum_number"] = len(NEG_GC_df)
    scripts.global_report.report_dict["pos_gc_insilico_spectrum_number"] = len(POS_GC_In_Silico_df)
    scripts.global_report.report_dict["neg_gc_insilico_spectrum_number"] = len(NEG_GC_In_Silico_df)


def format_parameters():
    """
    Formats and presents the parameters dictionary into a structured string format for better readability.
    The function extracts input and output files, output directory configuration, and a variety of processing
    parameters, presenting them in a categorized and human-readable manner. It is specifically tailored for
    an input dictionary expected to have a predefined structure for file types and processing configurations.

    Arguments:
        parameters_dict (dict): A dictionary containing input file paths, output configurations, and processing
            parameter settings. Each key represents a parameter, its value either being a specific configuration
            or control flag.

    Returns:
        str: A formatted string summarizing the provided configuration and parameters, including input files,
            output format, directory details, and parameter-specific configurations. Detailed sections include
            file classifications and processing behavior based on input.
    """
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
        
    OUTPUT_FORMAT:
        CSV: {"YES" if parameters_dict['csv'] else "NO"}
        MSP: {"YES" if parameters_dict['msp'] else "NO"}
        jSON: {"YES" if parameters_dict['json'] else "NO"}

    ===================== PARAMETERS =====================
    normalize_intensity: {"ON" if parameters_dict["normalize_intensity"] else "OFF"}
    remove_peak_above_precursormz: {"ON" if parameters_dict['remove_peak_above_precursormz'] else "OFF"}
    check_minimum_peak_requiered: {"ON" if parameters_dict['check_minimum_peak_requiered'] else "OFF"}
        n_peaks: {parameters_dict['check_minimum_peak_requiered_n_peaks']}
    reduce_peak_list: {"ON" if parameters_dict['reduce_peak_list'] else "OFF"}
        max_peaks: {parameters_dict['reduce_peak_list_max_peaks']}
    remove_spectrum_under_entropy_score: {"ON" if parameters_dict['remove_spectrum_under_entropy_score'] else "OFF"}
        entropy_score_value: {parameters_dict['remove_spectrum_under_entropy_score_value']}
    keep_mz_in_range: {"ON" if parameters_dict['keep_mz_in_range'] else "OFF"}
        from_mz: {parameters_dict['keep_mz_in_range_from_mz']}
        to_mz: {parameters_dict['keep_mz_in_range_to_mz']}
    check_minimum_of_high_peaks_requiered: {"ON" if parameters_dict['check_minimum_of_high_peaks_requiered'] else "OFF"}
        intensity_percent: {parameters_dict['check_minimum_of_high_peaks_requiered_intensity_percent']}
        no_peaks: {parameters_dict['check_minimum_of_high_peaks_requiered_no_peaks']}
    
    reset_updates: {"YES" if parameters_dict['reset_updates'] else "NO"}
    
    """

    return parameters_string


def format_fitered_out():
    """
    Generate a summary report of items that have been filtered out based on specific
    criteria. This function formats the filtering summary into a predefined structured
    string.

    Returns
    -------
    str
        A formatted string summarizing the items that were filtered out.

    """
    filtered_out_string = f""" 
    ======================= FILTERED OUT =======================
    No peaks list: {scripts.deletion_report.no_peaks_list}
    No smiles, no inchi, no inchikey: {scripts.deletion_report.no_smiles_no_inchi_no_inchikey}
    No precursor mz: {scripts.deletion_report.no_precursor_mz}
    No or bad adduct: {scripts.deletion_report.no_or_bad_adduct}
    Low entropy score: {scripts.deletion_report.low_entropy_score}
    Minimum peaks not required: {scripts.deletion_report.minimum_peaks_not_requiered}
    All peaks above precursor mz: {scripts.deletion_report.all_peaks_above_precursor_mz}
    No peaks in mz range: {scripts.deletion_report.no_peaks_in_mz_range}
    Minimum high peaks not required: {scripts.deletion_report.minimum_high_peaks_not_requiered}
    
    """

    return filtered_out_string


def format_spectrum_numbers():
    """
    Generate a summary string containing detailed spectrum information alongside
    a total spectrum count. Computes the total spectrum number by summing individual
    components and formats it into a structured multi-line string for better readability.

    Returns
    -------
    str
        A formatted string containing spectrum information including individual
        spectrum numbers categorized by type and mode (Positive LC Exp, Negative GC
        Insilico, etc.) and the total spectrum number.
    """
    total_spectrum_number = (
            scripts.global_report.report_dict["pos_lc_exp_spectrum_number"] +
            scripts.global_report.report_dict["neg_lc_exp_spectrum_number"] +
            scripts.global_report.report_dict["pos_lc_insilico_spectrum_number"] +
            scripts.global_report.report_dict["neg_lc_insilico_spectrum_number"] +
            scripts.global_report.report_dict["pos_gc_exp_spectrum_number"] +
            scripts.global_report.report_dict["neg_gc_exp_spectrum_number"] +
            scripts.global_report.report_dict["pos_gc_insilico_spectrum_number"] +
            scripts.global_report.report_dict["neg_gc_insilico_spectrum_number"]
    )

    spectrum_numbers_string = f""" 
    ================== SPECTRUM NUMBER ==================
    POS LC Exp: {scripts.global_report.report_dict["pos_lc_exp_spectrum_number"]}
    NEG LC Exp: {scripts.global_report.report_dict["neg_lc_exp_spectrum_number"]}
    POS LC InSilico: {scripts.global_report.report_dict["pos_lc_insilico_spectrum_number"]}
    NEG LC InSilico: {scripts.global_report.report_dict["neg_lc_insilico_spectrum_number"]}
    POS GC Exp: {scripts.global_report.report_dict["pos_gc_exp_spectrum_number"]}
    NEG GC Exp: {scripts.global_report.report_dict["neg_gc_exp_spectrum_number"]}
    POS GC InSilico: {scripts.global_report.report_dict["pos_gc_insilico_spectrum_number"]}
    NEG GC InSilico: {scripts.global_report.report_dict["neg_gc_insilico_spectrum_number"]}

    Total: {total_spectrum_number}
    """

    return spectrum_numbers_string


def format_unique_inchikeys():
    """
    Generates a detailed formatted string containing counts of unique InChIKeys from a global report
    dictionary. The unique InChIKeys are categorized based on spectrum type (LC/GC), ionization mode
    (POS/NEG), and source type (Exp/InSilico).

    Returns
    -------
    str
        A formatted multiline string containing unique InChIKey counts for each category and the
        total count.

    """
    unique_inchikeys_string = f""" 
    ================= UNIQUE INCHIKEYS ==================
    POS LC Exp: {scripts.global_report.report_dict["pos_lc_exp_spectrum_unique_inchikey"]}
    NEG LC Exp: {scripts.global_report.report_dict["neg_lc_exp_spectrum_unique_inchikey"]}
    POS LC InSilico: {scripts.global_report.report_dict["pos_lc_insilico_spectrum_unique_inchikey"]}
    NEG LC InSilico: {scripts.global_report.report_dict["neg_lc_insilico_spectrum_unique_inchikey"]}
    POS GC Exp: {scripts.global_report.report_dict["pos_gc_exp_spectrum_unique_inchikey"]}
    NEG GC Exp: {scripts.global_report.report_dict["neg_gc_exp_spectrum_unique_inchikey"]}
    POS GC InSilico: {scripts.global_report.report_dict["pos_gc_insilico_spectrum_unique_inchikey"]}
    NEG GC InSilico: {scripts.global_report.report_dict["neg_gc_insilico_spectrum_unique_inchikey"]}

    TOTAL Unique InChIKeys: {scripts.global_report.report_dict["TOTAL_unique_inchikey"]}
    
    """

    return unique_inchikeys_string


def format_report():
    """
    Generates a formatted report string by concatenating strings representing different
    sections of the report.

    This function combines information about parameters, filtered-out data, spectrum
    numbers, and unique InChIKeys into a single formatted report string.

    Returns:
        str: The concatenated formatted string containing the report details.

    Raises:
        None
    """
    parameters_string = format_parameters()
    filtered_out_string = format_fitered_out()
    spectrum_numbers_string = format_spectrum_numbers()
    unique_inchikeys_string = format_unique_inchikeys()

    return parameters_string + filtered_out_string + spectrum_numbers_string + unique_inchikeys_string


def report(output_directory, POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df):
    """
    Generate a formatted report based on input data and write it to a specified output directory.

    This function calculates spectra and unique InChIKeys from the given input lists.
    It then generates a formatted report and writes the content to a file named "report.txt"
    in the specified output directory.

    Parameters:
        output_directory (str): The directory where the output file, "report.txt", will be written.
        POS_LC (list): Positive liquid chromatography input data.
        POS_LC_insilico (list): In-silico generated data for positive liquid chromatography.
        POS_GC (list): Positive gas chromatography input data.
        POS_GC_insilico (list): In-silico generated data for positive gas chromatography.
        NEG_LC (list): Negative liquid chromatography input data.
        NEG_LC_insilico (list): In-silico generated data for negative liquid chromatography.
        NEG_GC (list): Negative gas chromatography input data.
        NEG_GC_insilico (list): In-silico generated data for negative gas chromatography.

    Raises:
        FileNotFoundError: If the specified output directory does not exist.
        IOError: If there is an error writing to the file in the specified output directory.

    Returns:
        None
    """
    # Calculs des spectres et des InChIKeys uniques
    calculate_spectrum_number(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df)
    calculate_unique_inchikeys(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df)

    # Création du rapport formaté
    formated_report = format_report()

    # Générer la date et l'heure au format désiré : DD_MM_YYYY__HH_MM_SS
    current_datetime = datetime.now().strftime("%d_%m_%Y__%H_%M_%S")

    # Définir le chemin complet pour le fichier avec date et heure
    report_file_path = os.path.join(output_directory, f"report_{current_datetime}.txt")

    # Écriture du contenu dans le fichier report.txt
    with open(report_file_path, "w") as report_file:
        report_file.write(formated_report)
