from scripts.GUI.utils.global_vars import parameters_dict
from scripts.normalizer.mols_calculation import *
from scripts.complete_from_pubchem_datas import *
from scripts.convertors.parsing_to_dict import *
from scripts.ontologies_completion import *
from scripts.convertors.csv_to_msp import *
from scripts.spectrum_normalizer import *
from scripts.duplicatas_remover import *
from scripts.splash_generator import *
from scripts.set_projects import *
import scripts.deletion_report
from scripts.splitter import *
from scripts.writers import *
from scripts.update import *
from scripts.report import *
import traceback
import time
import sys
import os


ordered_columns = ["FILENAME",
                   "PREDICTED",
                   "SPLASH",
                   "SPECTRUMID",
                   "RESOLUTION",
                   "SYNON",
                   "IONIZATION",
                   "MSLEVEL",
                   "FRAGMENTATIONMODE",
                   "NAME",
                   "PRECURSORMZ",
                   "EXACTMASS",
                   "AVERAGEMASS",
                   "PRECURSORTYPE",
                   "INSTRUMENTTYPE",
                   "INSTRUMENT",
                   "SMILES",
                   "INCHI",
                   "INCHIKEY",
                   "COLLISIONENERGY",
                   "FORMULA",
                   "RT",
                   "IONMODE",
                   "COMMENT",
                   "ENTROPY",
                   "CLASSYFIRE_SUPERCLASS",
                   "CLASSYFIRE_CLASS",
                   "CLASSYFIRE_SUBCLASS",
                   "NPCLASS_PATHWAY",
                   "NPCLASS_SUPERCLASS",
                   "NPCLASS_CLASS",
                   "NUM PEAKS",
                   "PEAKS_LIST"]

def check_interrupt(stop_flag):
    """Helper function to check the stop flag."""
    if stop_flag and stop_flag():
        print("Interruption requested by the user.")
        raise InterruptedError("Execution stopped by user request via stop_flag.")

def MAIN(progress_callback=None, total_items_callback=None, prefix_callback=None, item_type_callback=None, step_callback=None, completion_callback=None, deletion_callback=None, stop_flag=None):
    """
    Executes the FragHub multi-step processing pipeline, including file parsing, data transformations, and cleaning,
    to process and analyze spectral data.

    This function consolidates several processing stages, such as file conversion, cleaning, metadata completion,
    and data segmentation, to prepare spectral datasets for advanced analysis. It handles various file formats, removes
    duplicate records, ensures data consistency, and retrieves missing metadata when necessary. The processed data is
    further split into subsets based on predefined characteristics for downstream workflows.

    Parameters:
        progress_callback (Callable, optional): A callback function to monitor the progress of operations.
        total_items_callback (Callable, optional): A callback function that provides the total number of items to process.
        prefix_callback (Callable, optional): A callback function for updating prefix names during processing.
        item_type_callback (Callable, optional): A callback function to keep track of the type of items being processed.
        step_callback (Callable, optional): A callback function signaling the start or end of specific steps.
        completion_callback (Callable, optional): A callback function to notify about the process's overall completion.
        deletion_callback (Callable, optional): A callback function to report deletions during processing.
        stop_flag (multiprocessing.Event, optional): Used to check if the process should be interrupted.

    Returns:
        int: A return code where 0 indicates successful completion or early exit if no files are found to process.
    """
    output_directory = parameters_dict["output_directory"]

    try:
        check_interrupt(stop_flag)  # Verify before starting

        if parameters_dict['reset_updates'] == 1.0:
            reset_updates(output_directory)
            check_interrupt(stop_flag)

        init_project(output_directory)
        check_interrupt(stop_flag)

        start_time = time.time()

        input_path = parameters_dict["input_directory"]

        # STEP 1: Convert files to JSON if needed (Multithreaded)
        FINAL_MSP, FINAL_CSV, FINAL_JSON, FINAL_MGF = parsing_to_dict(input_path, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback, step_callback=step_callback)
        check_interrupt(stop_flag)

        files_to_process = False

        # Check if there are files to process
        if FINAL_MSP or FINAL_CSV or FINAL_JSON or FINAL_MGF:
            files_to_process = True

        # If no files are available to process, stop execution
        if not files_to_process:
            deletion_callback("-- THERE ARE NO FILES TO PROCESS. EXITING PROCESS --")
            time.sleep(0.01)
            if completion_callback:
                completion_callback("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
            return 0

        # STEP 2: Generate SPLASH KEY
        time.sleep(0.01)
        if step_callback:
            step_callback("-- GENERATING SPLASH UNIQUE ID --")
        time.sleep(0.01)
        FINAL_MSP, FINAL_CSV, FINAL_JSON, FINAL_MGF = generate_splash_id(FINAL_MSP, FINAL_CSV, FINAL_JSON, FINAL_MGF, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        check_interrupt(stop_flag)

        spectrum_list = []
        spectrum_list.extend(FINAL_MSP)
        del FINAL_MSP
        spectrum_list.extend(FINAL_CSV)
        del FINAL_CSV
        spectrum_list.extend(FINAL_JSON)
        del FINAL_JSON
        spectrum_list.extend(FINAL_MGF)
        del FINAL_MGF
        check_interrupt(stop_flag)

        # STEP 3: Remove duplicates
        spectrum_list = pd.DataFrame(spectrum_list)[ordered_columns]
        # Convert all columns to str except 'PEAKS_LIST'
        spectrum_list = spectrum_list.astype({col: str for col in ordered_columns if col != 'PEAKS_LIST'})
        time.sleep(0.01)
        if step_callback:
            step_callback("-- REMOVING DUPLICATES --")
        time.sleep(0.01)
        spectrum_list = remove_duplicatas(spectrum_list, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        deletion_callback(f"Duplicates removed: {scripts.deletion_report.duplicatas_removed}")
        check_interrupt(stop_flag)

        first_run = False
        update = False

        # STEP 4: Clean spectra (Multithreaded)
        time.sleep(0.01)
        if step_callback:
            step_callback("-- CHECKING FOR UPDATES --")
        time.sleep(0.01)
        spectrum_list, update_temp, first_run_temp = check_for_update_processing(spectrum_list, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
        deletion_callback(f"Previously cleaned: {scripts.deletion_report.previously_cleaned}")
        check_interrupt(stop_flag)

        if spectrum_list:

            if update_temp:
                update = True
            if first_run_temp:
                first_run = True
            time.sleep(0.01)
            if step_callback:
                step_callback("-- CLEANING SPECTRA --")
            time.sleep(0.01)
            spectrum_list = spectrum_cleaning_processing(spectrum_list, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
            deletion_callback(
                f"""
                No peaks list: {scripts.deletion_report.no_peaks_list}
                No SMILES, no InChI, no InChIKey: {scripts.deletion_report.no_smiles_no_inchi_no_inchikey}
                No precursor m/z: {scripts.deletion_report.no_precursor_mz}
                No or bad adduct: {scripts.deletion_report.no_or_bad_adduct}
                Low entropy score: {scripts.deletion_report.low_entropy_score}
                Minimum peaks not required: {scripts.deletion_report.minimum_peaks_not_requiered}
                All peaks above precursor m/z: {scripts.deletion_report.all_peaks_above_precursor_mz}
                No peaks in m/z range: {scripts.deletion_report.no_peaks_in_mz_range}
                Minimum high peaks not reached: {scripts.deletion_report.minimum_high_peaks_not_requiered}
                """
            )
            check_interrupt(stop_flag)

            if not spectrum_list:
                deletion_callback("-- THERE ARE NO SPECTRA TO PROCESS AFTER CLEANING. EXITING PROCESS --")
                time.sleep(0.01)
                if completion_callback:
                    completion_callback("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
                return 0

            spectrum_list = pd.DataFrame(spectrum_list)[ordered_columns].astype(str)

            # STEP 5: mols derivations and calculations
            time.sleep(0.01)
            if step_callback:
                step_callback("--  MOLS DERIVATION AND MASS CALCULATION --")
            time.sleep(0.01)
            spectrum_list = mols_derivation_and_calculation(spectrum_list, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
            deletion_callback(f"No smiles, no inchi, no inchikey (updated): {scripts.deletion_report.no_smiles_no_inchi_no_inchikey}")
            check_interrupt(stop_flag)

            # STEP 6: completing missing metadata from pubchem datas
            time.sleep(0.01)
            if step_callback:
                step_callback("--  COMPLETING FROM PUBCHEM DATAS --")
            time.sleep(0.01)
            spectrum_list = complete_from_pubchem_datas(spectrum_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)

            # STEP 7: completing missing names
            time.sleep(0.01)
            if step_callback:
                step_callback("--  ONTOLOGIES COMPLETION --")
            time.sleep(0.01)
            spectrum_list = ontologies_completion(spectrum_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
            check_interrupt(stop_flag)

            # STEP 8: SPLITTING
            # -- SPLITTING [POS / NEG] --
            time.sleep(0.01)
            if step_callback:
                step_callback("--  SPLITTING [POS / NEG] --")
            time.sleep(0.01)
            POS_df, NEG_df = split_pos_neg(spectrum_list, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
            check_interrupt(stop_flag)

            # -- SPLITTING [LC / GC] --
            time.sleep(0.01)
            if step_callback:
                step_callback("--  SPLITTING [LC / GC] --")
            time.sleep(0.01)
            POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df = split_LC_GC(POS_df, NEG_df, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
            check_interrupt(stop_flag)

            del POS_df
            del NEG_df

            # -- SPLITTING [EXP / In-Silico] --
            time.sleep(0.01)
            if step_callback:
                step_callback("--  SPLITTING [EXP / In-Silico] --")
            time.sleep(0.01)
            POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df = exp_in_silico_splitter(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
            check_interrupt(stop_flag)
            if parameters_dict["msp"] == 1.0:
                time.sleep(0.01)
                if step_callback:
                    step_callback("--  CONVERTING CSV TO MSP --")
                time.sleep(0.01)
                POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico = csv_to_msp(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
                check_interrupt(stop_flag)
            # STEP 9: writting output files
            if parameters_dict["csv"] == 1.0:
                time.sleep(0.01)
                if step_callback:
                    step_callback("--  WRITING CSV --")
                time.sleep(0.01)
                writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_In_Silico_df, POS_GC_In_Silico_df, NEG_LC_In_Silico_df, NEG_GC_In_Silico_df, first_run, output_directory, update, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
                check_interrupt(stop_flag)
            if parameters_dict["msp"] == 1.0:
                time.sleep(0.01)
                if step_callback:
                    step_callback("--  WRITING MSP --")
                time.sleep(0.01)
                writting_msp(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico, output_directory, update, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
                check_interrupt(stop_flag)
            if parameters_dict["json"] == 1.0:
                time.sleep(0.01)
                if step_callback:
                    step_callback("--  WRITING JSON --")
                time.sleep(0.01)
                writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_In_Silico_df, POS_GC_In_Silico_df, NEG_LC_In_Silico_df, NEG_GC_In_Silico_df, output_directory, progress_callback=progress_callback, total_items_callback=total_items_callback, prefix_callback=prefix_callback, item_type_callback=item_type_callback)
                check_interrupt(stop_flag)
            deletion_callback(
                f"Total deletions: {sum([scripts.deletion_report.duplicatas_removed, scripts.deletion_report.previously_cleaned, scripts.deletion_report.no_peaks_list, scripts.deletion_report.no_smiles_no_inchi_no_inchikey, scripts.deletion_report.no_precursor_mz, scripts.deletion_report.low_entropy_score, scripts.deletion_report.minimum_peaks_not_requiered, scripts.deletion_report.all_peaks_above_precursor_mz, scripts.deletion_report.no_peaks_in_mz_range, scripts.deletion_report.minimum_high_peaks_not_requiered])}"
            )
        else:
            deletion_callback("There is no new spectrums to process. Exiting code !")

        report(output_directory, POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df)
        check_interrupt(stop_flag)

        time.sleep(0.01)
        if completion_callback:
            completion_callback("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))

    except InterruptedError as ie:
        # Specifically handled in run_main_function, but we can log it here
        print(f"MAIN function interrupted: {ie}")

    except Exception as e:
        # In case of an unexpected error in MAIN, propagate it to run_main_function
        print(f"!!! Unhandled exception occurred in MAIN: {e}")
        traceback.print_exc()  # Display in the thread worker console
        raise  # Let run_main_function handle and emit the appropriate signal