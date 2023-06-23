from bs4 import BeautifulSoup as bs
import concurrent.futures
from tqdm import tqdm
import pandas as pd
import json
import lxml
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
        results = list(tqdm(executor.map(json_to_msp, FINAL_JSON), total=len(FINAL_JSON)))

    final = [res for res in results if res is not None]

    return final # returns the list of different worker executions.

def complete_HMDB(specrta_dict):
    HMDB_df = pd.read_csv("../datas/HMDB.csv",sep=";",encoding="UTF-8")
    HMDB_ID = specrta_dict["database-id"] # hmdb id retrieval for csv matching

    line = HMDB_df.loc[HMDB_df['EXTERNAL_ID'] == HMDB_ID]

    specrta_dict["inchikey"] = line["INCHIKEY"].values[0]
    specrta_dict["smiles"] = line["SMILES"].values[0]
    specrta_dict["inchi"] = line["INCHI"].values[0]
    specrta_dict["formula"] = line["MOLECULAR_FORMULA"].values[0]
    specrta_dict["PRECURSORMZ"] = line["ACC_MASS"].values[0]

    return specrta_dict

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

    soup = bs(xml_content, "lxml")

    # Starting retrieve all information contained in xml markdown
    # specrta_dict["sample-concentration"] = soup.find("sample-concentration")
    # if specrta_dict["sample-concentration"] != None:
    #     specrta_dict["sample-concentration"] = specrta_dict["sample-concentration"].contents
    #
    # specrta_dict["solvent"] = soup.find("solvent")
    # if specrta_dict["solvent"] != None:
    #     specrta_dict["solvent"] = specrta_dict["solvent"].contents
    #
    # specrta_dict["sample-mass"] = soup.find("sample-mass")
    # if specrta_dict["sample-mass"] != None:
    #     specrta_dict["sample-mass"] = specrta_dict["sample-mass"].contents
    #
    # specrta_dict["sample-assessment"] = soup.find("sample-assessment")
    # if specrta_dict["sample-assessment"] != None:
    #     specrta_dict["sample-assessment"] = specrta_dict["sample-assessment"].contents
    #
    # specrta_dict["spectra-assessment"] = soup.find("spectra-assessment")
    # if specrta_dict["spectra-assessment"] != None:
    #     specrta_dict["spectra-assessment"] = specrta_dict["spectra-assessment"].contents
    #
    # specrta_dict["sample-source"] = soup.find("sample-source")
    # if specrta_dict["sample-source"] != None:
    #     specrta_dict["sample-source"] = specrta_dict["sample-source"].contents
    #
    # specrta_dict["collection-date"] = soup.find("collection-date")
    # if specrta_dict["collection-date"] != None:
    #     specrta_dict["collection-date"] = specrta_dict["collection-date"].contents

    specrta_dict["filename"] = soup.find("filename")
    if specrta_dict["filename"] != None:
        specrta_dict["filename"] = specrta_dict["filename"].contents


    specrta_dict["instrument-type"] = soup.find("instrument-type")
    if specrta_dict["instrument-type"] != None:
        specrta_dict["instrument-type"] = specrta_dict["instrument-type"].contents

    # specrta_dict["peak-counter"] = soup.find("peak-counter")
    # if specrta_dict["peak-counter"] != None:
    #     specrta_dict["peak-counter"] = specrta_dict["peak-counter"].contents
    #
    # specrta_dict["created-at"] = soup.find("created-at")
    # if specrta_dict["created-at"] != None:
    #     specrta_dict["created-at"] = specrta_dict["created-at"].contents
    #
    # specrta_dict["updated-at"] = soup.find("updated-at")
    # if specrta_dict["updated-at"] != None:
    #     specrta_dict["updated-at"] = specrta_dict["updated-at"].contents

    # specrta_dict["mono-mass"] = soup.find("mono-mass")
    # if specrta_dict["mono-mass"] != None:
    #     specrta_dict["mono-mass"] = specrta_dict["mono-mass"].contents
    #
    # specrta_dict["energy-field"] = soup.find("energy-field")
    # if specrta_dict["energy-field"] != None:
    #     specrta_dict["energy-field"] = specrta_dict["energy-field"].contents
    #
    # specrta_dict["collision-energy-level"] = soup.find("collision-energy-level")
    # if specrta_dict["collision-energy-level"] != None:
    #     specrta_dict["collision-energy-level"] = specrta_dict["collision-energy-level"].contents

    specrta_dict["collision-energy-voltage"] = soup.find("collision-energy-voltage")
    if specrta_dict["collision-energy-voltage"] != None:
        specrta_dict["collision-energy-voltage"] = specrta_dict["collision-energy-voltage"].contents

    specrta_dict["ionization-mode"] = soup.find("ionization-mode")
    if specrta_dict["ionization-mode"] != None:
        specrta_dict["ionization-mode"] = specrta_dict["ionization-mode"].contents

    # specrta_dict["sample-concentration-units"] = soup.find("sample-concentration-units")
    # if specrta_dict["sample-concentration-units"] != None:
    #     specrta_dict["sample-concentration-units"] = specrta_dict["sample-concentration-units"].contents
    #
    # specrta_dict["sample-mass-units"] = soup.find("sample-mass-units")
    # if specrta_dict["sample-mass-units"] != None:
    #     specrta_dict["sample-mass-units"] = specrta_dict["sample-mass-units"].contents

    specrta_dict["predicted"] = soup.find("predicted")
    if specrta_dict["predicted"] != None:
        specrta_dict["predicted"] = specrta_dict["predicted"].contents

    # specrta_dict["structure-id"] = soup.find("structure-id")
    # if specrta_dict["structure-id"] != None:
    #     specrta_dict["structure-id"] = specrta_dict["structure-id"].contents
    #
    # specrta_dict["splash-key"] = soup.find("splash-key")
    # if specrta_dict["splash-key"] != None:
    #     specrta_dict["splash-key"] = specrta_dict["splash-key"].contents
    #
    # specrta_dict["chromatography-type"] = soup.find("chromatography-type")
    # if specrta_dict["chromatography-type"] != None:
    #     specrta_dict["chromatography-type"] = specrta_dict["chromatography-type"].contents
    #
    # specrta_dict["analyzer-type"] = soup.find("analyzer-type")
    # if specrta_dict["analyzer-type"] != None:
    #     specrta_dict["analyzer-type"] = specrta_dict["analyzer-type"].contents
    #
    # specrta_dict["ionization-type"] = soup.find("ionization-type")
    # if specrta_dict["ionization-type"] != None:
    #     specrta_dict["ionization-type"] = specrta_dict["ionization-type"].contents
    #
    # specrta_dict["charge-type"] = soup.find("charge-type")
    # if specrta_dict["charge-type"] != None:
    #     specrta_dict["charge-type"] = specrta_dict["charge-type"].contents
    #
    # specrta_dict["data-source"] = soup.find("data-source")
    # if specrta_dict["data-source"] != None:
    #     specrta_dict["data-source"] = specrta_dict["data-source"].contents
    #
    # specrta_dict["data-source-id"] = soup.find("data-source-id")
    # if specrta_dict["data-source-id"] != None:
    #     specrta_dict["data-source-id"] = specrta_dict["data-source-id"].contents

    specrta_dict["adduct"] = soup.find("adduct")
    if specrta_dict["adduct"] != None:
        specrta_dict["adduct"] = specrta_dict["adduct"].contents

    # specrta_dict["adduct-type"] = soup.find("adduct-type")
    # if specrta_dict["adduct-type"] != None:
    #     specrta_dict["adduct-type"] = specrta_dict["adduct-type"].contents
    #
    # specrta_dict["adduct-mass"] = soup.find("adduct-mass")
    # if specrta_dict["adduct-mass"] != None:
    #     specrta_dict["adduct-mass"] = specrta_dict["adduct-mass"].contents

    specrta_dict["database-id"] = soup.find("database-id")
    if specrta_dict["database-id"] != None:
        specrta_dict["database-id"] = specrta_dict["database-id"].contents

    # Correcting 0 charge
    if specrta_dict["ionization-mode"][0] == "Positive":
        specrta_dict["charge"] = "1"
    elif specrta_dict["ionization-mode"][0] == "Negative":
        specrta_dict["charge"] = "-1"

    peak_list = [[mz.contents for mz in soup.find_all("mass-charge")],
                 [intensity.contents for intensity in soup.find_all("intensity")]]

    specrta_dict_final = {}

    for key, value in specrta_dict.items():
        specrta_dict_final[key] = None
        if specrta_dict[key] != [] and specrta_dict[key] != None:
            specrta_dict_final[key] = specrta_dict[key][0]

    # Starting to write information from xml to msp format
    specrta_dict_final["inchikey"] = ""
    specrta_dict_final["smiles"] = ""
    specrta_dict_final["inchi"] = ""
    specrta_dict_final["formula"] = ""
    specrta_dict_final["PRECURSORMZ"] = ""
    specrta_dict_final["peak_list"] = ""

    for mass_charge, intensity in zip(peak_list[0], peak_list[1]):
        specrta_dict_final["peak_list"] = specrta_dict_final["peak_list"] + mass_charge[0] + " " + intensity[0] + "\n"

    # complete HMDB spectrum with PKON23 datas
    specrta_dict_final = complete_HMDB(specrta_dict_final)

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
        results = list(tqdm(executor.map(xml_to_msp, FINAL_XML), total=len(FINAL_XML)))

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
        print("CONVERTING JSON TO MSP")
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
        print("CONVERTING XML TO MSP")
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

            for spectrum in spectrum_list:
                if spectrum != "\n":
                    fields = re.findall(r"(.+?):(.*)\n", spectrum)
                    if first == True:
                        for element in fields:
                            dictionary[element[0]] = []
                        first = False

                    if first == False:
                        for element in fields:
                            dictionary[element[0]].append(element[1])

                    dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2).split("\n"))

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

            for spectrum in spectrum_list:
                if spectrum != "\n":
                    fields = re.findall(r"(.+?):(.*)", spectrum)
                    if first == True:
                        for element in fields:
                            dictionary[element[0]] = []
                        first = False

                    if first == False:
                        for element in fields:
                            dictionary[element[0]].append(element[1])

                    dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectrum).group(2).split("\n"))

    NEG_df = pd.DataFrame.from_dict(dictionary)
    NEG_df.to_csv("../OUTPUT/CSV/FINAL_NEG/NEG_clean.csv", sep=";", encoding="UTF-8", index=False)



