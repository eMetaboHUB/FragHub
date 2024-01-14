import json
import decimal
import ijson
import itertools

def decimal_default(obj):
    if isinstance(obj, decimal.Decimal):
        return float(obj)
    raise TypeError

with open('C:/Users/Axel/PycharmProjects/msp_v3/INPUT/CONVERTED/JSON_converted.json', 'r', encoding="UTF-8") as fichier:
    # Creation de l'iterator
    objects = ijson.items(fichier, 'item')

    # Obtenir le 100ème élement
    centieme_element = next(itertools.islice(objects, 99999, None), None)

    # Affichage du 100ème element en pretty JSON
    if centieme_element:
        print(json.dumps(centieme_element, default=decimal_default, indent=4))

# import ijson
# import json
# from decimal import Decimal
#
# # Cette fonction sera utilisée pour convertir les objets Decimal en float
# def decimal_default(obj):
#     if isinstance(obj, Decimal):
#         return float(obj)
#     raise TypeError
#
# def process_large_json(input_filename, output_filename, num_records):
#     with open(input_filename, 'r', encoding='utf-8') as input_file:
#         objects = ijson.items(input_file, 'item')
#
#         output_list = []
#
#         for idx, obj in enumerate(objects):
#             output_list.append(obj)
#             if idx + 1 == num_records:
#                 break
#
#         with open(output_filename, 'w') as output_file:
#             json.dump(output_list, output_file, default=decimal_default, indent=4)  # Utilisez la fonction decimal_default ici
#             print(f'{num_records} spectres ont été écrits dans {output_filename}')
#
#
# # Appeler la fonction
# process_large_json(r'C:\Users\Axel\PycharmProjects\msp_v3\INPUT\JSON\LITTLE_MoNA-export-Experimental_Spectra.json', r'C:\Users\Axel\PycharmProjects\msp_v3\INPUT\JSON\LITTLE_MoNA-export-Experimental_Spectra.json', 10000)