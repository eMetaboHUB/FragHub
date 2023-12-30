
from duplicatas_remover import *
from msp_normalizer import *
from set_parameters import *
from fraghubid_generator import *
import numpy as np

from rdkit import RDLogger
from msp_utilities import *
from converters import *
from splitter import *
from writers import *
import logging
import time
import sys
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages

if __name__ == "__main__":

    # Execution du GUI
    build_window()

    start_time = time.time()

    input_path = r"..\INPUT"
    output_path = r"..\OUTPUT"

    # STEP 1: convert files to msp if needed (Multithreaded)
    FINAL_JSON,FINAL_XML,FINAL_CSV = convert_to_msp(input_path)

    if FINAL_JSON != [] :
        with open(os.path.join("../INPUT/MSP/JSON_converted"+".msp"), "a", encoding="UTF-8") as temp:
            temp.write("\n\n".join(FINAL_JSON))

    if FINAL_XML != []:
        with open(os.path.join("../INPUT/MSP/XML_converted"+".msp"), "a", encoding="UTF-8") as temp:
            temp.write("\n\n".join(FINAL_XML))

    if FINAL_CSV != []:
        compteur = 1
        for elements in FINAL_CSV:
            with open(os.path.join(f"../INPUT/MSP/CSV_converted_{compteur}"+".msp"), "w", encoding="UTF-8") as temp:
                temp.write(elements)
            compteur += 1

    del FINAL_JSON
    del FINAL_XML
    del FINAL_CSV

    # STEP 2: split spectrums into a list AND matchms clean
    msp_dir = os.path.join(input_path, "MSP")

    # Check if there is msp file to process
    msp_to_process = False
    for files in os.listdir(msp_dir):
        if files.endswith(".msp"):
            msp_to_process = True
            break

    # If there is no msp to process: stop python execution
    if msp_to_process == False:
        sys.exit("There is no msp file to process into \"./INPUT/MSP\". Exiting code !")

    # generating FRAGHUBID
    print("-- GENERATING FragHub UNIQUE ID --")
    time.sleep(0.01)
    generate_fraghub_id(r"../INPUT/MSP")

    # Creating list of spectrums
    for files in os.listdir(msp_dir):
        if files.endswith(".msp"):
            msp_path = os.path.join(msp_dir, files)

            # correct_uncomplete_charge(msp_path)

            # STEP 3: Execute matchms (Multithreaded)
            print(f"-- CLEANING: {files} --")
            spectrum_list = load_spectrum_list(msp_path)
            spectrum_list = msp_cleaning_processing(spectrum_list)

            del spectrum_list

            # Write matchms clean msp into new msp file
            with open(os.path.join(r"../OUTPUT/MSP", file_name + "_clean" + ".msp"), "w", encoding="UTF-8") as clean:
                clean.write("\n\n".join(results))

            del results

    # # STEP 4: (All msp files were cleaned)
    clean_msp_path = os.path.join(output_path,"MSP")
    CONCATENATE_LIST = concatenate_clean_msp(clean_msp_path)  # ICI mettre un df en sortie pour démarrer le traitement csv

    CONCATENATE_DF = msp_to_csv(CONCATENATE_LIST)

    del CONCATENATE_LIST

    print("-- MOLS HARMONIZATION --")
    time.sleep(0.01)
    CONCATENATE_DF = mols_derivator(CONCATENATE_DF)  # REMOVE NO SMILES/INCHI désormais inclut dans cette fonction

    print("-- MASS CALCULATION --")
    time.sleep(0.01)
    CONCATENATE_DF = mass_calculator(CONCATENATE_DF)

    print("-- NAMES COMPLETION --")
    time.sleep(0.01)
    CONCATENATE_DF = names_completion(CONCATENATE_DF)

    print("-- SPLITTING [POS / NEG] --")
    time.sleep(0.01)
    POS, NEG = split_pos_neg(CONCATENATE_DF)

    # STEP 5: Split LC / GC
    time.sleep(0.01)
    print("-- SPLITTING [LC / GC] --")
    time.sleep(0.01)
    POS_LC,POS_GC,NEG_LC,NEG_GC = split_LC_GC(POS,NEG)

    del POS
    del NEG

    # STEP 5: EXP / In-Silico
    time.sleep(0.01)
    print("-- SPLITTING EXP / In-Silico --")
    time.sleep(0.01)
    POS_LC,POS_LC_In_Silico,POS_GC,POS_GC_In_Silico,NEG_LC,NEG_LC_In_Silico,NEG_GC,NEG_GC_In_Silico = exp_in_silico_splitter(POS_LC,POS_GC,NEG_LC,NEG_GC)


    # STEP 6: Remove duplicates spectrum when same peak_list for the same inchikey.
    print("-- REMOVING DUPLICATAS --")
    time.sleep(0.01)
    POS_LC,POS_LC_df,POS_LC_df_insilico,POS_LC_In_Silico,POS_GC,POS_GC_df,POS_GC_df_insilico,POS_GC_In_Silico,NEG_LC,NEG_LC_df,NEG_LC_df_insilico,NEG_LC_In_Silico,NEG_GC,NEG_GC_df,NEG_GC_df_insilico,NEG_GC_In_Silico = remove_duplicatas(POS_LC,POS_LC_In_Silico,POS_GC,POS_GC_In_Silico,NEG_LC,NEG_LC_In_Silico,NEG_GC,NEG_GC_In_Silico)

    print("-- WRITING MSP --")
    writting_msp(clean_msp_path, POS_LC, POS_GC, NEG_LC, NEG_GC, POS_LC_In_Silico, POS_GC_In_Silico, NEG_LC_In_Silico, NEG_GC_In_Silico)

    print("-- WRITING CSV --")
    writting_csv(POS_LC_df, POS_GC_df, NEG_LC_df, NEG_GC_df, POS_LC_df_insilico, POS_GC_df_insilico, NEG_LC_df_insilico, NEG_GC_df_insilico)

    print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))