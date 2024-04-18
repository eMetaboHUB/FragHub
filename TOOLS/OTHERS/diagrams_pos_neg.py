import os
import pandas as pd
import matplotlib.pyplot as plt



# AFTER
df_1 = pd.read_csv(r"C:\Users\Axel\Documents\PYTHON\MSP_V3\final_msp\All_DB_NEG_Jul2023.csv",sep=";",encoding="UTF-8")
df_2 = pd.read_csv(r"C:\Users\Axel\Documents\PYTHON\MSP_V3\final_msp\All_DB_POS_Jul2023.csv",sep=";",encoding="UTF-8")

df = pd.concat([df_1, df_2], ignore_index=True)

occurrences = df['charge'].value_counts()

# Créez un diagramme circulaire
plt.figure(figsize=(6, 6))
plt.pie(occurrences, labels=['POS', 'NEG'], autopct='%1.1f%%', startangle=140)
plt.axis('equal')  # Assurez-vous que le diagramme est un cercle
plt.title('Distribution of acquisition methods')

# Ajoutez une légende
plt.legend(loc='best', labels=['POS', 'NEG'])

# Affichez le diagramme circulaire
plt.show()

