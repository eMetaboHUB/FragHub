import re

texte = "MS2 Composite"

ms_level_pattern = re.compile("(?:ms)?(\d)", flags=re.IGNORECASE)

test = re.findall(ms_level_pattern, texte)

print(test)