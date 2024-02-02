import matplotlib
import matplotlib.cm as cm
import hashlib
import numpy as np

liste = ['Polyketides|Terpenoids', 'Polyketides', 'Alkaloids', 'Shikimates and Phenylpropanoids', 'Terpenoids',
         'Fatty acids', 'Amino acids and Peptides|Shikimates and Phenylpropanoids', 'Amino acids and Peptides',
         'Carbohydrates', 'Alkaloids|Shikimates and Phenylpropanoids', 'Carbohydrates|Shikimates and Phenylpropanoids',
         'Fatty acids|Polyketides', 'Polyketides|Shikimates and Phenylpropanoids',
         'Amino acids and Peptides|Fatty acids', 'Alkaloids|Amino acids and Peptides', 'Alkaloids|Terpenoids',
         'Amino acids and Peptides|Polyketides', 'Amino acids and Peptides|Fatty acids|Shikimates and Phenylpropanoids',
         'Alkaloids|Polyketides', 'Alkaloids|Amino acids and Peptides|Shikimates and Phenylpropanoids',
         'Shikimates and Phenylpropanoids|Terpenoids', 'Fatty acids|Terpenoids', 'Carbohydrates|Fatty acids',
         'Carbohydrates|Polyketides', 'Alkaloids|Carbohydrates', 'Amino acids and Peptides|Fatty acids|Polyketides',
         'Amino acids and Peptides|Carbohydrates', 'Alkaloids|Fatty acids', 'Fatty acids|Polyketides|Terpenoids']

# Récupérer valeurs uniques
valeurs_uniques = list(set(liste))

cmap = cm.get_cmap('nipy_spectral', len(valeurs_uniques))  # Colormap


def get_color(val):
    # Convertir la valeur en un nombre entre 0 et 1
    normalized = hashlib.md5(val.encode('utf-8')).hexdigest()
    normalized = int(normalized, 16) / 16 ** len(normalized)

    # Utiliser ce nombre pour obtenir une couleur
    return matplotlib.colors.rgb2hex(cmap(normalized))


# générer le dictionnaire couleur
couleurs = {val: get_color(val) for val in valeurs_uniques}

print(couleurs)