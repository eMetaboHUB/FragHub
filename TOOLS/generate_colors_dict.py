from matplotlib import colors
from colorsys import hsv_to_rgb

liste = ['Lipids and lipid-like molecules', 'Phenylpropanoids and polyketides', 'Organic acids and derivatives', 'Benzenoids', 'Organoheterocyclic compounds', 'Organohalogen compounds', 'Organic oxygen compounds', 'Organic nitrogen compounds', 'Nucleosides, nucleotides, and analogues', 'Lignans, neolignans and related compounds', 'Organosulfur compounds', 'Alkaloids and derivatives', 'Organophosphorus compounds', 'Hydrocarbon derivatives', 'Hydrocarbons', 'Homogeneous non-metal compounds', 'Organic 1,3-dipolar compounds', 'Organic Polymers', 'Organometallic compounds']
valeurs_uniques = list(set(liste))

# générer le dictionnaire couleur
couleurs = {val: colors.to_hex(hsv_to_rgb(i/len(valeurs_uniques), 1, 1)) for i, val in enumerate(valeurs_uniques)}

print(couleurs)