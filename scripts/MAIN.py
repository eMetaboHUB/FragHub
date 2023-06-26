from duplicatas_remover import *
from matchms_processing import *
from rdkit import RDLogger
from msp_utilities import *
from converters import *
from splitter import *
import logging
import sys
import os

RDLogger.DisableLog('rdApp.*') # Disable rdkit log (warning) messages
logger = logging.getLogger("matchms") # Disable matchms log messages
logger.disabled = True

import time

if __name__ == "__main__":
    start_time = time.time()
    
    input_path = r"..\INPUT"
    output_path = r"..\OUTPUT"

    # STEP 1: convert files to msp if needed (Multithreaded)
    FINAL_JSON,FINAL_XML = convert_to_msp(input_path)

    if FINAL_JSON != [] :
        with open(os.path.join("../INPUT/MSP/JSON_converted"+".msp"), "a", encoding="UTF-8") as temp:
            temp.write("\n\n".join(FINAL_JSON))

    if FINAL_XML != []:
        with open(os.path.join("../INPUT/MSP/XML_converted"+".msp"), "a", encoding="UTF-8") as temp:
            temp.write("\n\n".join(FINAL_XML))

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

    # Creating list of spectrums
    for files in os.listdir(msp_dir):
        if files.endswith(".msp"):
            msp_path = os.path.join(msp_dir, files)
            file_name = os.path.basename(msp_path.replace(".msp", ""))

            correct_uncomplete_charge(msp_path)

            # STEP 3: Execute matchms (Multithreaded)
            print("-- MATCHMS PROCESSING ON: ",file_name," --")
            spectrum_list = list(load_from_msp(msp_path))
            results = matchms_processing(spectrum_list,file_name)

            # Write matchms clean msp into new msp file
            with open(os.path.join(r"..\OUTPUT\CLEAN_MSP",file_name+"_clean"+".msp"), "w", encoding="UTF-8") as clean:
                clean.write("\n\n".join(results))

    # STEP 4: (All msp files were cleaned) --> Split POS and NEG spectrums
    clean_msp_path = os.path.join(output_path,"CLEAN_MSP")
    CONCATENATE_LIST = concatenate_clean_msp(clean_msp_path)

    print("-- SPLITTING POS / NEG --")
    time.sleep(0.01)
    POS, NEG = split_pos_neg(CONCATENATE_LIST)

    # STEP 5: Split LC / GC
    print("-- SPLITTING LC / GC --")
    time.sleep(0.01)
    POS_LC,POS_GC,NEG_LC,NEG_GC = split_LC_GC(POS,NEG)


    # STEP 6: Remove duplicates spectrum when same peak_list for the same inchikey.
    # print("-- REMOVING DUPLICATAS --")
    # time.sleep(0.01)
    # POS_LC_FILTERED,POS_GC_FILTERED,NEG_LC_FILTERED,NEG_GC_FILTERED = remove_duplicatas_public(POS_LC,POS_GC,NEG_LC,NEG_GC)

    print("REMOVING DUPLICATAS")
    time.sleep(0.01)
    POS_LC_FILTERED,POS_GC_FILTERED,NEG_LC_FILTERED,NEG_GC_FILTERED = remove_duplicatas_lrsv(POS_LC,POS_GC,NEG_LC,NEG_GC)

    POS_LC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(POS_LC_FILTERED))
    POS_GC_FULL = re.sub("\n{2,}", "\n\n\n", "\n\n".join(POS_GC_FILTERED))
    NEG_LC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG_LC_FILTERED))
    NEG_GC_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG_GC_FILTERED))

    with open(os.path.join(clean_msp_path,"FINAL_POS/POS_LC_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_LC_FULL)
    with open(os.path.join(clean_msp_path,"FINAL_POS/POS_GC_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_GC_FULL)
    with open(os.path.join(clean_msp_path, "FINAL_NEG/NEG_LC_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_LC_FULL)
    with open(os.path.join(clean_msp_path, "FINAL_NEG/NEG_GC_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_GC_FULL)

    print("-- CONVERTING MSP TO CSV --")
    msp_to_csv()

    print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))







