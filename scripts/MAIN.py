from msp_utilities import *
import os

input_path = r"..\INPUT"
output_path = r"..\OUTPUT"


# STEP 1: convert files to msp if needed
convert_to_msp(input_path)

# STEP 2: split spectrums into a list AND matchms clean
msp_dir = os.path.join(input_path,"MSP")
for files in os.listdir(msp_dir):
    if files.endswith(".msp"):
        msp_path = os.path.join(msp_dir,files)
        spectrums = split_spectrums(msp_path)

        # TO BE CONTINUED ...
