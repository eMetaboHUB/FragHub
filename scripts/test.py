from msp_parsers import *
from data_preparer import *
from set_parameters import *
import time
import os

# Execution de la fonction
build_window()

start_time = time.time()

# generating filepath list
spectrum_path_list = []
MSP_path = os.path.abspath(r"../INPUT/MSP")
for files in os.listdir(MSP_path):
    if files.endswith(".msp"):
        spectrum_path_list.append(os.path.join(MSP_path,files))

spectrum_path_list = [spectrums for spectrums in spectrum_path_list if spectrums != None]

for path in spectrum_path_list:
    spectrum_list = load_spectrum_list(path)
    break

# for spectrum in spectrum_list:
#     print("__",spectrum)
#     print("\n\n")

spectrum_list = msp_parsing_processing(spectrum_list)

df = pd.concat(spectrum_list, join='outer')

# df.to_csv(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test.csv",index=False, sep=";", quotechar='"', encoding="UTF-8")

df.to_excel(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test_msms_pos.xlsx",index=False)

# compteur = 1
# for spectrum in spectrum_list:
#     spectrum.to_excel(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\test_{compteur}.xlsx",index=False)
#     compteur += 1
#     break

print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))