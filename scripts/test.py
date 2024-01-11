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
    centieme_element = next(itertools.islice(objects, 9999, None), None)

    # Affichage du 100ème element en pretty JSON
    if centieme_element:
        print(json.dumps(centieme_element, default=decimal_default, indent=4))