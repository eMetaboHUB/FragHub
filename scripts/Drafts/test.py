from FragHub.scripts.msp_parsers import *
import time

start_time = time.time()

# generating filepath list
spectrum_path_list = []
MSP_path = os.path.abspath(r"../../INPUT/MSP")
for files in os.listdir(MSP_path):
    if files.endswith(".msp"):
        spectrum_path_list.append(os.path.join(MSP_path,files))

spectrum_path_list = [spectrums for spectrums in spectrum_path_list if spectrums != None]

for path in spectrum_path_list:
    spectrum_list = load_spectrum_list(path)
    print("spectrums loaded")
    spectrum_list = [spectrums for spectrums in spectrum_list if spectrums != None]
    spectrum_list = msp_parsing_processing(spectrum_list)

# print(spectrum_list)

compteur = 1
for spectrum in spectrum_list:
    spectrum.to_excel(rf"C:\Users\Axel\Documents\PYTHON\MSP_V3\FragHub\OUTPUT\MSP\TEST\test_{compteur}.xlsx",index=False)
    compteur += 1
    break

print("--- TOTAL TIME: %s ---" % time.strftime("%H:%M:%S", time.gmtime(time.time() - start_time)))