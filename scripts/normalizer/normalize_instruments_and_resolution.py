import json
import re


global instrument_tree
with open('../../datas/instruments_tree.json', 'r') as f:
    instrument_tree = json.load(f)

def clean_spectrum_instrument_info(metadata_dict):
    instrument_infos = metadata_dict["INSTRUMENT"] + " " + metadata_dict["INSTRUMENTTYPE"]
    instrument_infos = instrument_infos.lower()
    instrument_infos = re.sub("-tof", "tof", instrument_infos)
    instrument_infos = re.sub("q-", "q", instrument_infos)

def normalize_instruments_and_resolution(metadata_dict):
    """
    Normalize the instrument metadata in the given dictionary.

    :param metadata_dict: A dictionary containing the instrument metadata.
    :return: The normalized instrument metadata dictionary.
    """
    if metadata_dict["INSTRUMENT"]:
        metadata_dict_instrument = ' '.join(re.sub(r'[^\w\s]', ' ', metadata_dict["INSTRUMENT"].lower() + " " + metadata_dict["INSTRUMENTTYPE"].lower()).split())
        normalized_instrument_name = ''.join(sorted(list(metadata_dict_instrument)))

        if normalized_instrument_name in instruments_dict:
            metadata_dict["INSTRUMENT"] = f'{instruments_dict[normalized_instrument_name]["REF_INSTRUMENT"]}-{instruments_dict[normalized_instrument_name]["REF_MODELE"]}'
            metadata_dict["INSTRUMENTTYPE"] = f'{instruments_dict[normalized_instrument_name]["REF_SPECTRUM_TYPE"]}-{instruments_dict[normalized_instrument_name]["REF_IONISATION"]}-{instruments_dict[normalized_instrument_name]["REF_INSTRUMENT_TYPE"]}'
            metadata_dict["RESOLUTION"] = f'{instruments_dict[normalized_instrument_name]["REF_RESOLUTION"]}'
            return metadata_dict

    return metadata_dict