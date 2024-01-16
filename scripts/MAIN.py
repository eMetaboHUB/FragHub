from convertors.convert_to_json import *
from convertors.csv_to_msp import *
from fraghubid_generator import *
from duplicatas_remover import *
from name_completion import *
from msp_normalizer import *
from set_parameters import *
from rdkit import RDLogger
from splitter import *
from writers import *
from update import *
import json
import time
import sys
import os

RDLogger.DisableLog('rdApp.*')  # Disable rdkit log (warning) messages

ordered_columns = ["FILENAME",
                   "PREDICTED",
                   "FRAGHUBID",
                   "SPECTRUMID",
                   "RESOLUTION",
                   "SYNON",
                   "CHARGE",
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

    # Execution du GUI
    build_window()

    start_time = time.time()

    input_path = os.path.abspath(r"../INPUT")
    output_path = os.path.abspath(r"../OUTPUT")

    # STEP 1: convert files to json if needed (Multithreaded)
    FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON = convert_to_json(input_path)

    json_dir = os.path.join(input_path, "CONVERTED")

    files_to_process = False

    # Check if there is msp file to process
    if FINAL_MSP or FINAL_XML or FINAL_CSV or FINAL_JSON:
        files_to_process = True

    # If there is no msp to process: stop python execution
    if not files_to_process:
        sys.exit("There is no files to process. Exiting code !")

    # STEP 2: generating FRAGHUBID
    print("{:>70}".format("-- GENERATING FragHub UNIQUE ID --"))
    time.sleep(0.01)
    FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON = generate_fraghub_id(FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON)

    spectrum_list = []
    spectrum_list.extend(FINAL_MSP)
    del FINAL_MSP
    spectrum_list.extend(FINAL_XML)
    del FINAL_XML
    spectrum_list.extend(FINAL_CSV)
    del FINAL_CSV
    spectrum_list.extend(FINAL_JSON)
    del FINAL_JSON

    first_run = False
    update = False

    # STEP 3: cleaning spectrums (Multithreaded)
    time.sleep(0.01)
    print("{:>70}".format(f"-- CHECKING FOR UPDATES --"))
    spectrum_list, update_temp, first_run_temp = check_for_update_processing(spectrum_list)
    if update_temp:
        update = True
    if first_run_temp:
        first_run = True
    time.sleep(0.01)
    print("{:>70}".format(f"-- CLEANING SPECTRUMS --"))
    spectrum_list = spectrum_cleaning_processing(spectrum_list)

    spectrum_list = pd.DataFrame(spectrum_list)[ordered_columns].astype(str)

    # STEP 4: complete missing information into spectrum
    print("{:>70}".format("-- MOLS DERIVATION AND MASS CALCULATION --"))
    time.sleep(0.01)
    spectrum_list = mols_derivation_and_calculation(spectrum_list)

    print("{:>70}".format("-- NAMES COMPLETION --"))
    time.sleep(0.01)
    spectrum_list = names_completion(spectrum_list)

    # STEP 5: splitting POS/NEG -- LC/GC -- EXP/InSilico
    print("{:>70}".format("-- SPLITTING [POS / NEG] --"))
    time.sleep(0.01)
    POS_df, NEG_df = split_pos_neg(spectrum_list)

    time.sleep(0.01)
    print("{:>70}".format("-- SPLITTING [LC / GC] --"))
    time.sleep(0.01)
    POS_LC_df,POS_GC_df,NEG_LC_df,NEG_GC_df = split_LC_GC(POS_df, NEG_df)

    del POS_df
    del NEG_df

    time.sleep(0.01)
    print("{:>70}".format("-- SPLITTING EXP / In-Silico --"))
    time.sleep(0.01)
    POS_LC_df,POS_LC_In_Silico_df,POS_GC_df,POS_GC_In_Silico_df,NEG_LC_df,NEG_LC_In_Silico_df,NEG_GC_df,NEG_GC_In_Silico_df = exp_in_silico_splitter(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df)

    # STEP 6: Remove duplicates spectrum when same peak_list for the same inchikey.
    print("{:>70}".format("-- REMOVING DUPLICATAS --"))
    time.sleep(0.01)
    POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico = remove_duplicatas(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df, first_run, update)

    print("{:>70}".format("-- CONVERTING CSV TO MSP --"))
    time.sleep(0.01)
    POS_LC_df,POS_LC,POS_LC_df_insilico,POS_LC_insilico,POS_GC_df,POS_GC,POS_GC_df_insilico,POS_GC_insilico,NEG_LC_df,NEG_LC,NEG_LC_df_insilico,NEG_LC_insilico,NEG_GC_df,NEG_GC,NEG_GC_df_insilico,NEG_GC_insilico = csv_to_msp(POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico)

    # STEP 7: writting output files
    print("{:>70}".format("-- WRITING CSV --"))
    time.sleep(0.01)
    writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, first_run, update)

    print("{:>70}".format("-- WRITING MSP --"))
    time.sleep(0.01)
    writting_msp(POS_LC,POS_LC_insilico,POS_GC,POS_GC_insilico,NEG_LC,NEG_LC_insilico,NEG_GC,NEG_GC_insilico, update)

    print("{:>70}".format("-- WRITING JSON --"))
    time.sleep(0.01)
    writting_json(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico)

    time.sleep(0.01)
    print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))