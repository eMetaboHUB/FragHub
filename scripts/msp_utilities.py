from bs4 import BeautifulSoup as bs
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
                specrta_dict_final["peak_list"] = specrta_dict_final["peak_list"] + mass_charge[0] + " " + intensity[
                    0] + "\n"

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
