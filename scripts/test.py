import re

texte = "MS2-MS3 Composite"

ms_level_pattern = re.compile("(?:ms)?(\d)", flags=re.IGNORECASE)

ms_level = re.findall(ms_level_pattern, texte)
if ms_level:
    print(ms_level)
    if len(ms_level) == 1:
        test = ms_level[0]
    elif len(ms_level) >= 2:
        test = f"{ms_level[0]}-{ms_level[1]}"

print(test)