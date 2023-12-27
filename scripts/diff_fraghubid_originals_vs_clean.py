import pandas as pd
import os
import re


def diff_fraghubid_originals_vs_clean(fraghub_id_column):

    for files in os.listdir("../INPUT/MSP"):
        if files.endswith(".msp"):
            df = pd.DataFrame(columns=['FRAGHUBID'])  # define a DataFrame with a column 'FRAGHUBID'

            with open(os.path.abspath(os.path.join("../INPUT/MSP", files)), 'r', encoding='utf-8') as buffer:
                content = buffer.read()

            fragubid_list = re.findall("(?:FRAGHUBID: )(.*)", content)

            # Insert each FRAGHUBID into the DataFrame
            df['FRAGHUBID'] = fragubid_list
            # for fraghubid in fragubid_list:
            #     df = df.append(pd.DataFrame({'FRAGHUBID': [fraghubid]}), ignore_index=True)

            df_3 = df[~df['FRAGHUBID'].isin(fraghub_id_column['FRAGHUBID'])]

            df_3.to_csv(rf"C:\Users\Axel\PycharmProjects\msp_v3\OUTPUT\MSP\TEST\diff_{files}.csv", mode='w', sep=";", quotechar='"', encoding="UTF-8", index=False)
