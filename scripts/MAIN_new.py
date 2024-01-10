from convertors.convert_to_json import *
from fraghubid_generator import *
from duplicatas_remover import *
from name_completion import *
from msp_normalizer import *
from set_parameters import *
from rdkit import RDLogger
from converters import *
from splitter import *
from writers import *
from update import *
import json
import time
import sys
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

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

    input_path = r"../INPUT"
    output_path = r"../OUTPUT"

    # STEP 1: convert files to json if needed (Multithreaded)
    FINAL_MSP,FINAL_XML,FINAL_CSV = convert_to_json(input_path)

    if FINAL_MSP != []:
        with open(os.path.join("../INPUT/JSON/MSP_converted.json"), "w", encoding="UTF-8") as buffer:
            json.dump(FINAL_MSP, buffer, ensure_ascii=False, indent=4)

    if FINAL_XML != []:
        with open(os.path.join("../INPUT/JSON/XML_converted.json"), "w", encoding="UTF-8") as buffer:
            json.dump(FINAL_XML, buffer, ensure_ascii=False, indent=4)

    if FINAL_CSV != []:
        with open(os.path.join("../INPUT/JSON/CSV_converted.json"), "w", encoding="UTF-8") as buffer:
            json.dump(FINAL_CSV, buffer, ensure_ascii=False, indent=4)

    del FINAL_MSP
    del FINAL_XML
    del FINAL_CSV

    json_dir = os.path.join(input_path, "JSON")

    # Check if there is msp file to process
    for root, dirs, files in os.walk(json_dir):
        for file in files:
            if file.endswith(".json"):
                json_to_process = True
                break

    # If there is no msp to process: stop python execution
    if json_to_process == False:
        sys.exit("There is no json file to process into \"./INPUT/JSON\". Exiting code !")
    #
    # # STEP 2: generating FRAGHUBID
    # print("{:>80}".format("-- GENERATING FragHub UNIQUE ID --"))
    # time.sleep(0.01)
    # generate_fraghub_id(r"../INPUT/MSP")
    #
    # CONCATENATED_SPECTRUMS_RESULTS = []
    # first_run = False
    # update = False
    #
    # for root, dirs, files in os.walk(msp_dir):
    #     for file in files:
    #         if file.endswith(".msp"):
    #             msp_path = os.path.join(root, file)
    #
    #             # STEP 3: cleaning spectrums (Multithreaded)
    #             print("{:>80}".format(f"-- CLEANING: {file} --"))
    #             spectrum_list = load_spectrum_list(msp_path)
    #             final_spectrum_list, update_temp, first_run_temp = check_for_update_processing(spectrum_list)
    #             if update_temp:
    #                 update = True
    #             if first_run_temp:
    #                 first_run = True
    #             spectrum_list = msp_cleaning_processing(spectrum_list)
    #
    #             CONCATENATED_SPECTRUMS_RESULTS.extend(spectrum_list)
    #
    #             del spectrum_list
    #
    # CONCATENATED_SPECTRUMS_DATAFRAME = pd.DataFrame(CONCATENATED_SPECTRUMS_RESULTS)[ordered_columns].astype(str)
    #
    # del CONCATENATED_SPECTRUMS_RESULTS
    #
    # # STEP 4: complete missing information into spectrum
    # print("{:>80}".format("-- MOLS DERIVATION AND MASS CALCULATION --"))
    # time.sleep(0.01)
    # CONCATENATED_SPECTRUMS_DATAFRAME = mols_derivation_and_calculation(CONCATENATED_SPECTRUMS_DATAFRAME)
    #
    # print("{:>80}".format("-- NAMES COMPLETION --"))
    # time.sleep(0.01)
    # CONCATENATED_SPECTRUMS_DATAFRAME = names_completion(CONCATENATED_SPECTRUMS_DATAFRAME)
    #
    # # STEP 5: splitting POS/NEG -- LC/GC -- EXP/InSilico
    # print("{:>80}".format("-- SPLITTING [POS / NEG] --"))
    # time.sleep(0.01)
    # POS_df, NEG_df = split_pos_neg(CONCATENATED_SPECTRUMS_DATAFRAME)
    #
    # time.sleep(0.01)
    # print("{:>80}".format("-- SPLITTING [LC / GC] --"))
    # time.sleep(0.01)
    # POS_LC_df,POS_GC_df,NEG_LC_df,NEG_GC_df = split_LC_GC(POS_df, NEG_df)
    #
    # del POS_df
    # del NEG_df
    #
    # time.sleep(0.01)
    # print("{:>80}".format("-- SPLITTING EXP / In-Silico --"))
    # time.sleep(0.01)
    # POS_LC_df,POS_LC_In_Silico_df,POS_GC_df,POS_GC_In_Silico_df,NEG_LC_df,NEG_LC_In_Silico_df,NEG_GC_df,NEG_GC_In_Silico_df = exp_in_silico_splitter(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df)
    #
    # # STEP 6: Remove duplicates spectrum when same peak_list for the same inchikey.
    # print("{:>80}".format("-- REMOVING DUPLICATAS --"))
    # time.sleep(0.01)
    # POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico = remove_duplicatas(POS_LC_df, POS_LC_In_Silico_df, POS_GC_df, POS_GC_In_Silico_df, NEG_LC_df, NEG_LC_In_Silico_df, NEG_GC_df, NEG_GC_In_Silico_df, first_run, update)
    #
    # print("{:>80}".format("-- CONVERTING CSV TO MSP --"))
    # time.sleep(0.01)
    # POS_LC_df,POS_LC,POS_LC_df_insilico,POS_LC_insilico,POS_GC_df,POS_GC,POS_GC_df_insilico,POS_GC_insilico,NEG_LC_df,NEG_LC,NEG_LC_df_insilico,NEG_LC_insilico,NEG_GC_df,NEG_GC,NEG_GC_df_insilico,NEG_GC_insilico = csv_and_msp(POS_LC_df,POS_LC_df_insilico,POS_GC_df,POS_GC_df_insilico,NEG_LC_df,NEG_LC_df_insilico,NEG_GC_df,NEG_GC_df_insilico)
    #
    # # STEP 7: writting output files
    # print("{:>80}".format("-- WRITING CSV --"))
    # time.sleep(0.01)
    # writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico, first_run, update)
    #
    # print("{:>80}".format("-- WRITING MSP --"))
    # time.sleep(0.01)
    # writting_msp(POS_LC,POS_LC_insilico,POS_GC,POS_GC_insilico,NEG_LC,NEG_LC_insilico,NEG_GC,NEG_GC_insilico, update)

    time.sleep(0.01)
    print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))