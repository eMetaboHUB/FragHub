from normalizer.mols_calculation import *
from convertors.convert_to_json import *
from convertors.csv_to_msp import *
from fraghubid_generator import *
from duplicatas_remover import *
from name_completion import *
from msp_normalizer import *
from set_parameters import *
from splitter import *
from writers import *
from update import *
import time
import sys
import os

from TRACES import *

ordered_columns = ["FILENAME",
                   "PREDICTED",
                   "FRAGHUBID",
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
                   "RETENTIONTIME",
                   "IONMODE",
                   "COMMENT",
                   "NUM PEAKS",
                   "PEAKS_LIST"]

if __name__ == "__main__":

    # GUI execution
    build_window()

    profile_name = parameters_dict["selected_profile"]

    if parameters_dict['reset_updates'] == 1.0:
        reset_updates(profile_name)

    init_profile(profile_name)

    start_time = time.time()

    input_path = os.path.abspath(r"../INPUT")
    output_path = os.path.abspath(rf"../OUTPUT/{profile_name}")

    # STEP 1: convert files to json if needed (Multithreaded)
    FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF = convert_to_json(input_path)

    files_to_process = False

    # Check if there is msp file to process
    if FINAL_MSP or FINAL_XML or FINAL_CSV or FINAL_JSON or FINAL_MGF:
        files_to_process = True

    # If there is no msp to process: stop python execution
    if not files_to_process:
        sys.exit("There is no files to process. Exiting code !")

    # STEP 2: generating FRAGHUBID
    time.sleep(0.01)
    print("{:>70}".format("-- GENERATING FragHub UNIQUE ID --"))
    time.sleep(0.01)
    FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF = generate_fraghub_id(FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF)

    spectrum_list = []
    spectrum_list.extend(FINAL_MSP)
    del FINAL_MSP
    spectrum_list.extend(FINAL_XML)
    del FINAL_XML
    spectrum_list.extend(FINAL_CSV)
    del FINAL_CSV
    spectrum_list.extend(FINAL_JSON)
    del FINAL_JSON
    spectrum_list.extend(FINAL_MGF)
    del FINAL_MGF

    first_run = False
    update = False

    # STEP 3: cleaning spectrums (Multithreaded)
    time.sleep(0.01)
    print("{:>70}".format(f"-- CHECKING FOR UPDATES --"))
    time.sleep(0.01)
    spectrum_list, update_temp, first_run_temp = check_for_update_processing(spectrum_list, profile_name)

    if not spectrum_list:
        sys.exit("There is no new spectrums to clean from databases. Exiting code !")

    if update_temp:
        update = True
    if first_run_temp:
        first_run = True
    time.sleep(0.01)
    print("{:>70}".format(f"-- CLEANING SPECTRUMS --"))
    time.sleep(0.01)
    spectrum_list = spectrum_cleaning_processing(spectrum_list)

    spectrum_list_TRACES = spectrum_list

    spectrum_list = [spectrum[0] for spectrum in spectrum_list if spectrum[0] is not None]

    if not spectrum_list:
        sys.exit("There is no spectrums to process after cleaning. Exiting code !")

    spectrum_list = pd.DataFrame(spectrum_list)[ordered_columns].astype(str)

    # STEP 4: mols derivations and calculations
    time.sleep(0.01)
    print("{:>70}".format("-- MOLS DERIVATION AND MASS CALCULATION --"))
    time.sleep(0.01)
    spectrum_list = mols_derivation_and_calculation(spectrum_list)

    # STEP 5: completing missing names
    time.sleep(0.01)
    print("{:>70}".format("-- NAMES COMPLETION --"))
    time.sleep(0.01)
    spectrum_list = names_completion(spectrum_list)

    # STEP 6: SPLITTING
    # -- SPLITTING [POS / NEG] --
    time.sleep(0.01)
    print("{:>70}".format("-- SPLITTING [POS / NEG] --"))
    time.sleep(0.01)
    POS_df, NEG_df = split_pos_neg(spectrum_list)

    # -- SPLITTING [LC / GC] --
    time.sleep(0.01)
    print("{:>70}".format("-- SPLITTING [LC / GC] --"))
    time.sleep(0.01)
    POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df = split_LC_GC(POS_df, NEG_df)

    del POS_df
    del NEG_df

    # -- SPLITTING [EXP / In-Silico] --
    time.sleep(0.01)
    print("{:>70}".format("-- SPLITTING [EXP / In-Silico] --"))
    time.sleep(0.01)
    POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df = exp_in_silico_splitter(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df)

    # STEP 7: Remove duplicates spectrum when same peak_list for the same inchikey.
    time.sleep(0.01)
    print("{:>70}".format("-- REMOVING DUPLICATAS --"))
    time.sleep(0.01)
    POS_LC_df, POS_LC_df_insilico, POS_GC_df, POS_GC_df_insilico, NEG_LC_df, NEG_LC_df_insilico, NEG_GC_df, NEG_GC_df_insilico = remove_duplicatas(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df, first_run, profile_name, update)

    if parameters_dict["msp"] == 1.0:
        time.sleep(0.01)
        print("{:>70}".format("-- CONVERTING CSV TO MSP --"))
        time.sleep(0.01)
        POS_LC_df, POS_LC, POS_LC_df_insilico, POS_LC_insilico, POS_GC_df, POS_GC, POS_GC_df_insilico, POS_GC_insilico, NEG_LC_df, NEG_LC, NEG_LC_df_insilico, NEG_LC_insilico, NEG_GC_df, NEG_GC, NEG_GC_df_insilico, NEG_GC_insilico = csv_to_msp(POS_LC_df, POS_LC_df_insilico, POS_GC_df, POS_GC_df_insilico, NEG_LC_df, NEG_LC_df_insilico, NEG_GC_df, NEG_GC_df_insilico)

    # STEP 8: writting output files
    if parameters_dict["csv"] == 1.0:
        time.sleep(0.01)
        print("{:>70}".format("-- WRITING CSV --"))
        time.sleep(0.01)
        writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, first_run, profile_name, update)

    if parameters_dict["msp"] == 1.0:
        time.sleep(0.01)
        print("{:>70}".format("-- WRITING MSP --"))
        time.sleep(0.01)
        writting_msp(POS_LC, POS_LC_insilico, POS_GC, POS_GC_insilico, NEG_LC, NEG_LC_insilico, NEG_GC, NEG_GC_insilico, profile_name, update)

    if parameters_dict["json"] == 1.0:
        time.sleep(0.01)
        print("{:>70}".format("-- WRITING JSON --"))
        time.sleep(0.01)
        writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, profile_name)

    time.sleep(0.01)
    print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))
