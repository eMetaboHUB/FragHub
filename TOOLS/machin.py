import os
import re

with open(r"C:\Users\Axel\Documents\pharmakon_db\pkon23\sous_ensemble\test\genus(1).csv","r") as buffer:
    content = buffer.read()

content = re.sub("\"","",content)

with open(r"C:\Users\Axel\Documents\pharmakon_db\pkon23\sous_ensemble\test\genus_final.csv","w") as buffer:
    buffer.write(content)