from convertors.convert_to_json import *


global In_Silico_pattern
In_Silico_pattern = re.compile(r"in.silico|insilico|predicted|theoretical|Annotation.level.3", flags=re.IGNORECASE)



def in_filename(filename):
    """
    :param filename: The name of the file to be checked for the presence of the string "MSMS_Public" and matching a specific pattern.
    :return: Returns True if the filename does not contain "MSMS_Public" and matches a specific pattern defined by the In_Silico_pattern. Returns False otherwise.

    """
    if "MSMS_Public" not in filename:
        if re.search(In_Silico_pattern, filename):
            return True

    return False



def normalize_predicted(metadata_dict):
    """
    Normalize the predicted field in the given metadata dictionary.

    :param metadata_dict: A dictionary containing metadata information.
    :return: The updated metadata dictionary with the normalized predicted field.
    """
    comment_field = metadata_dict["COMMENT"]
    predicted = metadata_dict["PREDICTED"]
    filename = metadata_dict["FILENAME"]

    if re.search(In_Silico_pattern, comment_field) or predicted == "true" or in_filename(filename):
        metadata_dict["PREDICTED"] = "true"
        return metadata_dict
    else:
        metadata_dict["PREDICTED"] = "false"
        return metadata_dict



def convert_keys(dict_list):
    """
    Convert keys in metadata_dict based on the provided keys_dict and keys_list.

    :param metadata_dict: A dictionary containing metadata information.
    :return: A dictionary with converted keys based on the provided keys_dict and keys_list.
    """
    output = []
    for metadata_dict in dict_list:
        converted = {keys_dict[key.lower()]: val for key, val in metadata_dict.items() if key.lower() in keys_dict and keys_dict[key.lower()] in keys_list}

        converted.update({key: "" for key in keys_list if key not in converted})
        output.append(converted)

    return output



original_db_path = r"C:\Users\Axel\Documents\MSP_DB\ORIGINALS_msp_DB\DB Janvier 2024\DB_publi"

FINAL_MSP, FINAL_XML, FINAL_CSV, FINAL_JSON, FINAL_MGF = convert_to_json(original_db_path)

MSP_unique_keys = []
for sub_list in FINAL_MSP:
    MSP_unique_keys.extend(list(sub_list))

MSP_unique_keys = list(set(MSP_unique_keys))

XML_unique_keys = []
for sub_list in FINAL_XML:
    XML_unique_keys.extend(list(sub_list))

XML_unique_keys = list(set(XML_unique_keys))


CSV_unique_keys = []
for sub_list in FINAL_CSV:
    CSV_unique_keys.extend(list(sub_list))

CSV_unique_keys = list(set(CSV_unique_keys))


JSON_unique_keys = []
for sub_list in FINAL_JSON:
    JSON_unique_keys.extend(list(sub_list))

JSON_unique_keys = list(set(JSON_unique_keys))


MGF_unique_keys = []
for sub_list in FINAL_MGF:
    MGF_unique_keys.extend(list(sub_list))

MGF_unique_keys = list(set(MGF_unique_keys))

TOTAL_UNIQUE_KEYS = []
TOTAL_UNIQUE_KEYS.extend(MSP_unique_keys)
TOTAL_UNIQUE_KEYS.extend(XML_unique_keys)
TOTAL_UNIQUE_KEYS.extend(CSV_unique_keys)
TOTAL_UNIQUE_KEYS.extend(JSON_unique_keys)
TOTAL_UNIQUE_KEYS.extend(MGF_unique_keys)

print(TOTAL_UNIQUE_KEYS)

print(len(MSP_unique_keys))
print(len(XML_unique_keys))
print(len(CSV_unique_keys))
print(len(JSON_unique_keys))
print(len(MGF_unique_keys))

print(len(TOTAL_UNIQUE_KEYS))