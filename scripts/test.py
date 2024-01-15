# import json
# import decimal
# import ijson
# import itertools
#
# def decimal_default(obj):
#     if isinstance(obj, decimal.Decimal):
#         return float(obj)
#     raise TypeError
#
# with open('C:/Users/Axel/PycharmProjects/msp_v3/INPUT/CONVERTED/JSON_converted.json', 'r', encoding="UTF-8") as fichier:
#     # Creation de l'iterator
#     objects = ijson.items(fichier, 'item')
#
#     # Obtenir le 100ème élement
#     centieme_element = next(itertools.islice(objects, 99999, None), None)
#
#     # Affichage du 100ème element en pretty JSON
#     if centieme_element:
#         print(json.dumps(centieme_element, default=decimal_default, indent=4))
#
# # import ijson
# # import json
# # from decimal import Decimal
# #
# # # Cette fonction sera utilisée pour convertir les objets Decimal en float
# # def decimal_default(obj):
# #     if isinstance(obj, Decimal):
# #         return float(obj)
# #     raise TypeError
# #
# # def process_large_json(input_filename, output_filename, num_records):
# #     with open(input_filename, 'r', encoding='utf-8') as input_file:
# #         objects = ijson.items(input_file, 'item')
# #
# #         output_list = []
# #
# #         for idx, obj in enumerate(objects):
# #             output_list.append(obj)
# #             if idx + 1 == num_records:
# #                 break
# #
# #         with open(output_filename, 'w') as output_file:
# #             json.dump(output_list, output_file, default=decimal_default, indent=4)  # Utilisez la fonction decimal_default ici
# #             print(f'{num_records} spectres ont été écrits dans {output_filename}')
# #
# #
# # # Appeler la fonction
# # process_large_json(r'C:\Users\Axel\PycharmProjects\msp_v3\INPUT\JSON\LITTLE_MoNA-export-Experimental_Spectra.json', r'C:\Users\Axel\PycharmProjects\msp_v3\INPUT\JSON\LITTLE_MoNA-export-Experimental_Spectra.json', 10000)

peak_list = [[164.0571, 1.394777], [165.06506, 14.984178], [166.07202, 6.46953], [168.08833, 1.711479], [179.08055, 2.615611], [181.05972, 3.025871], [182.06694, 2.607008], [192.088, 4.237026], [195.07462, 1.073237], [196.08304, 11.603429], [208.08308, 32.228283], [209.08618, 2.857035], [218.06519, 54.508158], [219.06767, 3.9749], [220.08301, 5.686917], [234.03813, 3.629701], [238.0936, 56.634468], [239.09691, 5.741493], [256.10373, 82.664218], [257.10706, 8.149285], [278.08569, 100.0], [279.08871, 9.704564], [280.09116, 1.21465], [294.06, 8.641813]]

peaks = eval(peak_list)

print(peaks)