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

            matchms_treatment(spectrum_list)
            # TO BE CONTINUED ... to sth with clean_msp_list

            # Write matchms clean msp into new msp file
            with open(os.path.join(r"..\OUTPUT\CLEAN_MSP)",file_name+"_clean"+".msp"), "w", encoding="UTF-8") as clean:
                clean.write("\n\n".join(clean_msp_list))


