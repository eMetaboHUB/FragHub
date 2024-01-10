import re

whole_text = """
NAME: 3-Des-Microcystein_LR
PRECURSORMZ: 981.54
PRECURSORTYPE: M+H
FORMULA: C48H72N10O12
Ontology: 
INCHIKEY: IYDKWWDUBYWQGF-NNAZGLEUSA-N
INCHI: N/A
SMILES: CC(C)CC1NC(=O)C(C)NC(=O)C(=C)N(C)C(=O)CCC(NC(=O)C(C)C(NC(=O)C(CCCNC(N)=N)NC(=O)C(C)C(NC1=O)C(O)=O)\C=C\C(\C)=C\C(C)C(O)Cc1ccccc1)C(O)=O
RETENTIONTIME: CCS: 
IONMODE: Positive
INSTRUMENTTYPE: LC-ESI-qTof
INSTRUMENT: qTof
COLLISIONENERGY: 
Comment: DB#=CCMSLIB00000001547; origin=GNPS
Num Peaks: 218
289.286377    8068.0
295.545288    22507.0
298.489624    3925.0
317.324951    18742.0
319.655945    8604.0
"""

# L'expression régulière modifiée pour prendre en compte les espaces dans la clé et les caractères spéciaux comme '#'
pattern = re.compile(r"([^:\n]*?):\s*([^:\n]*)(?:\n|$)", re.MULTILINE)
matches = pattern.findall(whole_text)

for match in matches:
    print("Key: " + match[0].strip())  # Les fonctions strip() sont utilisées pour retirer les espaces en début et en fin de chaîne
    print("Value: " + match[1].strip())
    print("\n")