import random
import os


paths = "../../ORIGINALS_msp_DB/"

n = 100

for dir in os.listdir(paths):
    dir_path = os.path.join(paths,dir)
    if os.path.isdir(dir_path):
        for files in os.listdir(dir_path):
            msp_path = os.path.join(dir_path,files)
            if files.endswith(".msp"):
                with open(msp_path,"r",encoding="UTF-8") as msp_buffer:
                    spectrum_list = msp_buffer.read().split("\n\n")

                random_n_list = random.sample(spectrum_list, n)

                with open(os.path.join("../../SOUS_DB","SOUS_"+files.replace(".msp","")+".msp"),"w",encoding="UTF-8") as msp_buffer:
                    msp_buffer.write("\n\n".join(random_n_list))



