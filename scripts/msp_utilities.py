from bs4 import BeautifulSoup as bs
import concurrent.futures
import pandas as pd
import numpy as np
import json
import lxml
import re
import os


def json_to_msp(json_path):
    for files in os.listdir(json_path):
        if files.endswith(".json"):
            file_name = os.path.basename(os.path.join(json_path, files)).replace(".json", "")
            with  open(os.path.join(json_path, files), "r", encoding="UTF-8") as f:
                lines = f.readlines()

            data = [json.loads(line) for line in lines]  # returns JSON object as a list of dictionary

            SPECTRUM = ""  # Creating empty spectrum string

            for spectras in data:  # Creating spectrum to msp format from json format
                for key, value in spectras.items():
                    if key != "peaks":
                        SPECTRUM = SPECTRUM + key + ": " + str(value) + "\n"
                    else:
                        SPECTRUM = SPECTRUM + "num peaks: " + str(len(spectras["peaks"])) + "\n"
                        for fragments in spectras["peaks"]:
                            SPECTRUM = SPECTRUM + str(fragments[0]) + " " + str(fragments[1]) + "\n"
                SPECTRUM = SPECTRUM + "\n\n"

                with open(os.path.join("../INPUT/MSP",file_name+"converted"+".msp"), "w", encoding="UTF-8") as temp:  # Writing spectrum to msp format into msp directory
                    temp.write(SPECTRUM)

def xml_to_msp(xml_path):
    for files in os.listdir(xml_path):
        if files.endswith(".xml"):
            file_name = os.path.basename(os.path.join(xml_path, files).replace(".xml", ""))
            with open(os.path.join(xml_path, files), "r", encoding="UTF-8") as xml_file:
                xml_content = xml_file.read()
            specrta_dict = {}
            soup = bs(xml_content, "lxml")

            # Starting retrieve all information contained in xml markdown
            specrta_dict["sample-concentration"] = soup.find("sample-concentration")
            if specrta_dict["sample-concentration"] != None:
                specrta_dict["sample-concentration"] = specrta_dict["sample-concentration"].contents

            specrta_dict["solvent"] = soup.find("solvent")
            if specrta_dict["solvent"] != None:
                specrta_dict["solvent"] = specrta_dict["solvent"].contents

            specrta_dict["sample-mass"] = soup.find("sample-mass")
            if specrta_dict["sample-mass"] != None:
                specrta_dict["sample-mass"] = specrta_dict["sample-mass"].contents

            specrta_dict["sample-source"] = soup.find("sample-source")
            if specrta_dict["sample-source"] != None:
                specrta_dict["sample-source"] = specrta_dict["sample-source"].contents

            specrta_dict["instrument-type"] = soup.find("instrument-type")
            if specrta_dict["instrument-type"] != None:
                specrta_dict["instrument-type"] = specrta_dict["instrument-type"].contents

            specrta_dict["mono-mass"] = soup.find("mono-mass")
            if specrta_dict["mono-mass"] != None:
                specrta_dict["mono-mass"] = specrta_dict["mono-mass"].contents

            specrta_dict["energy-field"] = soup.find("energy-field")
            if specrta_dict["energy-field"] != None:
                specrta_dict["energy-field"] = specrta_dict["energy-field"].contents

            specrta_dict["collision-energy-level"] = soup.find("collision-energy-level")
            if specrta_dict["collision-energy-level"] != None:
                specrta_dict["collision-energy-level"] = specrta_dict["collision-energy-level"].contents

            specrta_dict["structure-id"] = soup.find("structure-id")
            if specrta_dict["structure-id"] != None:
                specrta_dict["structure-id"] = specrta_dict["structure-id"].contents

            specrta_dict["chromatography-type"] = soup.find("chromatography-type")
            if specrta_dict["chromatography-type"] != None:
                specrta_dict["chromatography-type"] = specrta_dict["chromatography-type"].contents

            specrta_dict["analyzer-type"] = soup.find("analyzer-type")
            if specrta_dict["analyzer-type"] != None:
                specrta_dict["analyzer-type"] = specrta_dict["analyzer-type"].contents

            specrta_dict["ionization-type"] = soup.find("ionization-type")
            if specrta_dict["ionization-type"] != None:
                specrta_dict["ionization-type"] = specrta_dict["ionization-type"].contents

            specrta_dict["data-source"] = soup.find("data-source")
            if specrta_dict["data-source"] != None:
                specrta_dict["data-source"] = specrta_dict["data-source"].contents

            specrta_dict["adduct"] = soup.find("adduct")
            if specrta_dict["adduct"] != None:
                specrta_dict["adduct"] = specrta_dict["adduct"].contents

            specrta_dict["adduct-type"] = soup.find("adduct-type")
            if specrta_dict["adduct-type"] != None:
                specrta_dict["adduct-type"] = specrta_dict["adduct-type"].contents

            specrta_dict["adduct-mass"] = soup.find("adduct-mass")
            if specrta_dict["adduct-mass"] != None:
                specrta_dict["adduct-mass"] = specrta_dict["adduct-mass"].contents

            specrta_dict["database-id"] = soup.find("database-id")
            if specrta_dict["database-id"] != None:
                specrta_dict["database-id"] = specrta_dict["database-id"].contents

            specrta_dict["database"] = soup.find("database")
            if specrta_dict["database"] != None:
                specrta_dict["database"] = specrta_dict["database"].contents

            peak_list = [[mz.contents for mz in soup.find_all("mass-charge")],
                         [intensity.contents for intensity in soup.find_all("intensity")]]

            specrta_dict_final = {}

            for key, value in specrta_dict.items():
                specrta_dict_final[key] = None
                if specrta_dict[key] != [] and specrta_dict[key] != None:
                    specrta_dict_final[key] = specrta_dict[key][0]

            # Starting to write information from xml to msp format
            specrta_dict_final["peak_list"] = ""

            for mass_charge, intensity in zip(peak_list[0], peak_list[1]):
                specrta_dict_final["peak_list"] = specrta_dict_final["peak_list"] + mass_charge[0] + " " + intensity[0] + "\n"

            MSP = ""
            for key, value in specrta_dict_final.items():
                if key == "peak_list":
                    MSP = MSP + "num peaks: " + str(len(peak_list[0])) + "\n"
                    MSP = MSP + str(specrta_dict_final[key]) + "\n"
                else:
                    MSP = MSP + key + ": " + str(specrta_dict_final[key]) + "\n"

            with open(os.path.join("../INPUT/MSP",file_name+"converted"+".msp"), "a", encoding="UTF-8") as temp:
                temp.write(MSP)
                temp.write("\n")


def convert_to_msp(input_path):
    # JSON
    json_to_do = False
    json_path = os.path.join(input_path,"JSON")
    # check if there is a json file into the directory
    for files in os.listdir(json_path):
        if files.endswith(".json"):
            json_to_do = True
            break
    if json_to_do == True:
        json_to_msp(json_path)

    # XML
    xml_to_do = False
    xml_path = os.path.join(input_path,"XML")
    # check if there is a xml file into the directory
    for files in os.listdir(xml_path):
        if files.endswith(".xml"):
            xml_to_do = True
            break
    if xml_to_do == True:
        xml_to_msp(xml_path)

def split_spectrums(msp_path):
    with open(msp_path,"r",encoding="UTF-8") as file_buffer:
        file_content = file_buffer.read()

    spectrums = file_content.split("\n\n")  # split spectrums into a list

    return spectrums

def concatenate_clean_msp(clean_msp_path):
    CONCATENATE_LIST = []
    # append all spectrums of all cleaned files into CONCATENATE_LIST
    for files in os.listdir(clean_msp_path):
        if files.endswith("_clean.msp"):
            with open(os.path.join(clean_msp_path,files),"r",encoding="UTF-8") as buffer:
                temp = buffer.read().split("\n\n")

            CONCATENATE_LIST.extend(temp)

    return CONCATENATE_LIST

def split_pos_neg(CONCATENATE_LIST):
    POS = []
    NEG = []
    for spectrum in CONCATENATE_LIST:
        SEARCH = re.search("(CHARGE): (.*)", spectrum)
        if SEARCH != None:
            if int(SEARCH.group(2)) > 0:
                POS.append(spectrum)
            else:
                NEG.append(spectrum)
    return POS, NEG


def harmonize_fields_names(file_path):
    df = pd.read_csv("./data/colomnsNames.csv", sep=";")
    dict = df.to_dict()

    # re-structure dict
    for column, elements in dict.items():
        dict[column] = [element for key, element in elements.items() if (element != column) and (str(element) != 'nan')]

    with open(file_path, "r", encoding="UTF-8") as file:  # open the temporary file
        file_content = file.read()

    # replace non_normalized fields by normalized fields
    for column, lists in dict.items():
        for non_normalize in lists:
            file_content = re.sub(non_normalize, column, file_content)

    with open(file_path, "r+", encoding="UTF-8") as file:  # open the temporary file
        file.truncate(0)

    return file_content

def msp2csv(clean_msp_path):
    # POS
    dictionary = {"file_name":[],"PEAKS_LIST":[]}

    POS_DIR = os.path.join(clean_msp_path,"FINAL_POS")

    for files in os.listdir(POS_DIR):
        if files.endswith(".msp"):
            with open(os.path.join(POS_DIR,files),"r",encoding="UTF-8") as msp_file_buffer:
                msp_file = msp_file_buffer.read()

            spectras_list = msp_file.split("\n\n") # une liste de spectres

            for spectras in spectras_list[:-1]:
                fields = re.findall(r"(.+?):(.*)",spectras)
                dico_temp = {}
                for element in fields:
                    dico_temp[element[0]] = element[1]

                dictionary["file_name"].append(files)

                for key,value in dico_temp.items():
                    if key not in dictionary.keys():
                        if len(dictionary["file_name"]) < 1:
                            dictionary[key] = []
                            dictionary[key].append(value)
                        else:
                            dictionary[key] = [np.nan for i in range(len(dictionary["file_name"])-1)]
                            dictionary[key].append(value)
                    else:
                        dictionary[key].append(value)

                # vérification inverse
                for key in dictionary.keys():
                    if key not in dico_temp.keys() and key != "file_name" and key != "PEAKS_LIST":
                        dictionary[key].append(np.nan)

                dictionary["PEAKS_LIST"].append(re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)",spectras).group(2).split("\n"))

    df_POS = pd.DataFrame.from_dict(dictionary)

    df_POS.to_csv(os.path.join(r"..\OUTPUT\CSV\FINAL_POS","POS.csv"),sep=";",encoding="UTF-8",index=False)

    # NEG
    dictionary = {"file_name": [], "PEAKS_LIST": []}

    NEG_DIR = os.path.join(clean_msp_path, "FINAL_NEG")

    for files in os.listdir(NEG_DIR):
        if files.endswith(".msp"):
            with open(os.path.join(NEG_DIR, files), "r", encoding="UTF-8") as msp_file_buffer:
                msp_file = msp_file_buffer.read()

            spectras_list = msp_file.split("\n\n")  # une liste de spectres

            for spectras in spectras_list[:-1]:
                fields = re.findall(r"(.+?):(.*)", spectras)
                dico_temp = {}
                for element in fields:
                    dico_temp[element[0]] = element[1]

                dictionary["file_name"].append(files)

                for key, value in dico_temp.items():
                    if key not in dictionary.keys():
                        if len(dictionary["file_name"]) < 1:
                            dictionary[key] = []
                            dictionary[key].append(value)
                        else:
                            dictionary[key] = [np.nan for i in range(len(dictionary["file_name"]) - 1)]
                            dictionary[key].append(value)
                    else:
                        dictionary[key].append(value)

                # vérification inverse
                for key in dictionary.keys():
                    if key not in dico_temp.keys() and key != "file_name" and key != "PEAKS_LIST":
                        dictionary[key].append(np.nan)

                dictionary["PEAKS_LIST"].append(
                    re.search("(NUM PEAKS: [0-9]*)\n([\s\S]*)", spectras).group(2).split("\n"))

    df_NEG = pd.DataFrame.from_dict(dictionary)

    df_NEG.to_csv(os.path.join(r"..\OUTPUT\CSV\FINAL_NEG", "NEG.csv"), sep=";", encoding="UTF-8", index=False)



