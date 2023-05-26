from tqdm.notebook import tqdm as tqdm
from msp_utilities import *
from matchms_treatment import *
import sys
import os

if __name__ == "__main__":
    input_path = r"..\INPUT"
    output_path = r"..\OUTPUT"

    # STEP 1: convert files to msp if needed
    convert_to_msp(input_path)

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
            spectrum_list = split_spectrums(msp_path)

            # STEP 3: Execute multithreaded matchms
            spectrum_list = tqdm(list(load_from_msp(msp_path)))
            matchms_treatment(spectrum_list)

            clean_content = harmonize_fields_names("./temp/temp.msp")

            # Write matchms clean msp into new msp file
            with open(os.path.join(r"..\OUTPUT\CLEAN_MSP",file_name+"_clean"+".msp"), "w", encoding="UTF-8") as clean:
                clean.write(clean_content)

    # STEP 4: (All msp files were cleaned) --> Split POS and NEG spectrums
    clean_msp_path = os.path.join(output_path,"CLEAN_MSP")
    CONCATENATE_LIST = concatenate_clean_msp(clean_msp_path)

    POS, NEG = split_pos_neg(CONCATENATE_LIST) # Multithreaded

    POS_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(POS))
    NEG_FULL = re.sub("\n{2,}","\n\n\n","\n\n".join(NEG))

    with open(os.path.join(clean_msp_path,"FINAL_POS/POS_clean.msp"),"w",encoding="UTF-8") as pos:
        pos.write(POS_FULL)
    with open(os.path.join(clean_msp_path, "FINAL_NEG/NEG_clean.msp"), "w", encoding="UTF-8") as neg:
        neg.write(NEG_FULL)

    msp2csv(clean_msp_path)






