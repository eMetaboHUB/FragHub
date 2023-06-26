import concurrent.futures
from tqdm import tqdm
import pandas as pd
import json
import time
import os
import re

def concatenate_json(json_path):
    JSON_LIST = []
    for files in os.listdir(json_path):
        if files.endswith(".json"):
            file_name = os.path.basename(os.path.join(json_path, files)).replace(".json", "")
            with  open(os.path.join(json_path, files), "r", encoding="UTF-8") as f:
                lines = f.readlines()

            data = [json.loads(line) for line in lines]  # returns JSON object as a list of dictionary
            # Add filename to json
            for dicts in data:
                dicts["filename"] = file_name

            JSON_LIST.extend(data)

    return JSON_LIST

def json_to_msp(json_spectrum):
    SPECTRUM = ""  # Creating empty spectrum string
    SPECTRUM = SPECTRUM + "FILENAME: " + json_spectrum["filename"] + "\n"
    for key, value in json_spectrum.items():
        if key != "peaks":
            SPECTRUM = SPECTRUM + key + ": " + str(value) + "\n"
        else:
            SPECTRUM = SPECTRUM + "num peaks: " + str(len(json_spectrum["peaks"])) + "\n"
            for fragments in json_spectrum["peaks"]:
                SPECTRUM = SPECTRUM + str(fragments[0]) + " " + str(fragments[1]) + "\n"

    return SPECTRUM

def JSON_convert_processing(FINAL_JSON):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(json_to_msp, FINAL_JSON), total=len(FINAL_JSON), unit="spectrums", colour="green"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.

def concatenate_xml(xml_path):
    FINAL_XML = []
    for files in os.listdir(xml_path):
        if files.endswith(".xml"):
            file_name = os.path.basename(os.path.join(xml_path, files).replace(".xml", ""))
            with open(os.path.join(xml_path, files), "r", encoding="UTF-8") as xml_file:
                xml_content = xml_file.read()
            # Add filename to xml
            xml_content = re.sub("</sample-mass>\n",f"</sample-mass>\n  <filename>{file_name}</filename>\n",xml_content)

            FINAL_XML.extend([xml_content])

    return FINAL_XML

def xml_to_msp(xml_content):
    specrta_dict = {}

    if re.search("<filename>(.*)</filename>", xml_content):
        specrta_dict["filename"] = re.search("<filename>(.*)</filename>", xml_content).group(1)

    if re.search("<inchikey>(.*)</inchikey>", xml_content):
        specrta_dict["inchikey"] = re.search("<inchikey>(.*)</inchikey>", xml_content).group(1)

    if re.search("<SMILES>(.*)</SMILES>", xml_content):
        specrta_dict["SMILES"] = re.search("<SMILES>(.*)</SMILES>", xml_content).group(1)

    if re.search("<INCHI>(.*)</INCHI>", xml_content):
        specrta_dict["INCHI"] = re.search("<INCHI>(.*)</INCHI>", xml_content).group(1)

    if re.search("<MOLECULAR_FORMULA>(.*)</MOLECULAR_FORMULA>", xml_content):
        specrta_dict["MOLECULAR_FORMULA"] = re.search("<MOLECULAR_FORMULA>(.*)</MOLECULAR_FORMULA>", xml_content).group(1)

    if re.search("<ACC_MASS>(.*)</ACC_MASS>", xml_content):
        specrta_dict["PRECURSORMZ"] = re.search("<ACC_MASS>(.*)</ACC_MASS>", xml_content).group(1)

    if re.search("<instrument-type>(.*)</instrument-type>", xml_content):
        specrta_dict["instrument-type"] = re.search("<instrument-type>(.*)</instrument-type>", xml_content).group(1)

    if re.search("<collision-energy-voltage>(.*)</collision-energy-voltage>", xml_content):
        specrta_dict["collision-energy-voltage"] = re.search("<collision-energy-voltage>(.*)</collision-energy-voltage>", xml_content).group(1)

    if re.search("<ionization-mode>(.*)</ionization-mode>", xml_content):
        specrta_dict["ionization-mode"] = re.search("<ionization-mode>(.*)</ionization-mode>", xml_content).group(1)

    if re.search("<predicted>(.*)</predicted>", xml_content):
        specrta_dict["predicted"] = re.search("<predicted>(.*)</predicted>", xml_content).group(1)

    if re.search("<database-id>(.*)</database-id>", xml_content):
        specrta_dict["database-id"] = re.search("<database-id>(.*)</database-id>", xml_content).group(1)

    # Correcting 0 charge
    if specrta_dict["ionization-mode"][0] == "Positive":
        specrta_dict["charge"] = "1"
    elif specrta_dict["ionization-mode"][0] == "Negative":
        specrta_dict["charge"] = "-1"

    peak_list = [re.findall("<mass-charge>(.*)</mass-charge>",xml_content),re.findall("<intensity>(.*)</intensity>",xml_content)]

    specrta_dict_final = {}

    for key, value in specrta_dict.items():
        specrta_dict_final[key] = None
        if specrta_dict[key] != [] and specrta_dict[key] != None:
            specrta_dict_final[key] = specrta_dict[key]

    # Starting to write information from xml to msp format
    specrta_dict_final["peak_list"] = ""

    for mass_charge, intensity in zip(peak_list[0], peak_list[1]):
        specrta_dict_final["peak_list"] = specrta_dict_final["peak_list"] + mass_charge + " " + intensity + "\n"

    SPECTRUM = ""

    for key, value in specrta_dict_final.items():
        if key == "peak_list":
            SPECTRUM = SPECTRUM + "NUM PEAKS: " + str(len(peak_list[0])) + "\n"
            SPECTRUM = SPECTRUM + str(specrta_dict_final[key])
        else:
            SPECTRUM = SPECTRUM + key + ": " + str(specrta_dict_final[key]) + "\n"

    return SPECTRUM

def XML_convert_processing(FINAL_XML):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        results = list(tqdm(executor.map(xml_to_msp, FINAL_XML), total=len(FINAL_XML), unit="spectrums", colour="green"))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.


def convert_to_msp(input_path):
    # JSON
    FINAL_JSON = []
    json_to_do = False
    json_path = os.path.join(input_path,"JSON")
    # check if there is a json file into the directory
    for files in os.listdir(json_path):
        if files.endswith(".json"):
            json_to_do = True
    if json_to_do == True:
        print("-- CONVERTING JSON TO MSP --")
        # Concatenate all JSON to a list
        FINAL_JSON = concatenate_json(json_path)
        # Convert all JSON spectrum to MSP spectrum (Multithreaded)
        FINAL_JSON = JSON_convert_processing(FINAL_JSON)


    # XML
    FINAL_XML = []
    xml_to_do = False
    xml_path = os.path.join(input_path,"XML")
    # check if there is a xml file into the directory
    for files in os.listdir(xml_path):
        if files.endswith(".xml"):
            xml_to_do = True
    if xml_to_do == True:
        print("-- CONVERTING XML TO MSP --")
        # Concatenate all XML to a list
        FINAL_XML = concatenate_xml(xml_path)
        # Convert all XML spectrum to MSP spectrum (Multithreaded)
        FINAL_XML = XML_convert_processing(FINAL_XML)

    return FINAL_JSON,FINAL_XML

def msp_to_csv(clean_msp_path):
    # POS
    POS_DIR = os.path.join(clean_msp_path, "FINAL_POS")

    dictionary = {"PEAKS_LIST": []}

    first = True

    for files in os.listdir(POS_DIR):
        if files.endswith(".msp"):
            with open(os.path.join(POS_DIR, files), "r", encoding="UTF-8") as msp_file_buffer:
                msp_file = msp_file_buffer.read()

            spectrum_list = msp_file.split("\n\n")  # une liste de spectres

            empty = False
            if len(spectrum_list) == 0:
                empty = True

            print("POS")
            time.sleep(0.01)
            for spectrum in list(tqdm(spectrum_list, total=len(spectrum_list), unit="spectrums", colour="green")):
                if spectrum != "\n":
                    fields = re.findall(r"(.+?):(.*)\n", spectrum)
                    if first == True:
                        for element in fields:
                            dictionary[element[0]] = []
                        first = False

                    if first == False:
                        for element in fields:
                            dictionary[element[0]].append(element[1])

                    if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectrum):
                        dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2).split("\n"))

    if empty == False:
        POS_df = pd.DataFrame.from_dict(dictionary)
        POS_df.to_csv("../OUTPUT/CSV/FINAL_POS/POS_clean.csv", sep=";", encoding="UTF-8", index=False)

    # NEG
    NEG_DIR = os.path.join(clean_msp_path, "FINAL_NEG")

    dictionary = {"PEAKS_LIST": []}

    first = True

    for files in os.listdir(NEG_DIR):
        if files.endswith(".msp"):
            with open(os.path.join(NEG_DIR, files), "r", encoding="UTF-8") as msp_file_buffer:
                msp_file = msp_file_buffer.read()

            spectrum_list = msp_file.split("\n\n")  # une liste de spectres

            empty = False
            if len(spectrum_list) == 0:
                empty = True

            print("NEG")
            time.sleep(0.01)
            for spectrum in list(tqdm(spectrum_list, total=len(spectrum_list), unit="spectrums", colour="green")):
                if spectrum != "\n":
                    fields = re.findall(r"(.+?):(.*)", spectrum)
                    if first == True:
                        for element in fields:
                            dictionary[element[0]] = []
                        first = False

                    if first == False:
                        for element in fields:
                            dictionary[element[0]].append(element[1])

                    if re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum):
                        dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2).split("\n"))

    if empty == False:
        NEG_df = pd.DataFrame.from_dict(dictionary)
        NEG_df.to_csv("../OUTPUT/CSV/FINAL_NEG/NEG_clean.csv", sep=";", encoding="UTF-8", index=False)



