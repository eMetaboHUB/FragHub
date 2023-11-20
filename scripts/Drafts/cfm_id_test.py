from multiprocessing import Process, Manager
import pandas as pd
import subprocess
import os
import re


def cfm_id_treatment(mols,final_list,total,files):
    metadata = re.search("([\s\S]*(.*\n)*NUM PEAKS: (.*\n))", mols).group(1) # Recover matadata

    SMILES = None
    INCHI = None
    CHARGE = None

    if re.search("SMILES: (.*)\n", mols): # Searching for SMILES in the spectrum metadata
        SMILES = re.search("SMILES: (.*)\n", mols).group(1)
    elif re.search("INCHI: (.*)\n", mols): # Searching for INCHI in the spectrum metadata
        INCHI = re.search("INCHI: (.*)\n", mols).group(1)

    if re.search("CHARGE: (.*)\n", mols): # Searching for CHARGE in the spectrum metadata
        CHARGE = re.search("CHARGE: (.*)\n", mols).group(1)

    completed = None
    if SMILES != None: # do cmf-id with the smiles
        if CHARGE == "1": # POS
            cmd = "docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cd /cfmid/public/; cfm-predict \'" + SMILES + "\' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 stdout\""
            completed = subprocess.run(["powershell", "-Command", cmd], capture_output=True).stdout.decode("UTF-8") # Running cfm-id docker into container via Windows Powershell
        elif CHARGE == "-1": # NEG
            cmd = "docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cd /cfmid/public/; cfm-predict \'" + SMILES + "\' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 stdout\""
            completed = subprocess.run(["powershell", "-Command", cmd], capture_output=True).stdout.decode("UTF-8") # Running cfm-id docker into container via Windows Powershell
    elif INCHI != None: # do cfm-id with the inchi
        if CHARGE == "1": # POS
            cmd = "docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cd /cfmid/public/; cfm-predict \'" + INCHI + "\' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 stdout\""
            completed = subprocess.run(["powershell", "-Command", cmd], capture_output=True).stdout.decode("UTF-8") # Running cfm-id docker into container via Windows Powershell
        elif CHARGE == "-1": # NEG
            cmd = "docker run --rm=true -v ${pwd}:/cfmid/public/ -i wishartlab/cfmid:latest sh -c \"cd /cfmid/public/; cfm-predict \'" + INCHI + "\' 0.001 /trained_models_cfmid4.0/[M-H]-/param_output.log /trained_models_cfmid4.0/[M-H]-/param_config.txt 1 stdout\""
            completed = subprocess.run(["powershell", "-Command", cmd], capture_output=True).stdout.decode("UTF-8") # Running cfm-id docker into container via Windows Powershell
    else: # NO smiles or inchi for this spectrum metadata
        print("NO smiles or inchi for this molécule.")
        pass

    if completed != None:

        # Moving Some metadata fields into a "COMMENT" field.
        temp = [None,None]

        if re.search("(#In-silico .*)", completed):
            temp[0] = re.search("(#In-silico .*)", completed).group(1)
        if re.search("(#PREDICTED BY .*)", completed):
            temp[1] = re.search("(#PREDICTED BY .*)", completed).group(1)

        comment = None

        if temp != [None,None]:
            comment = "COMMENT: " + temp[0] + " " + temp[1] + "\n"
        # Moving DONE

        if comment != None:
            metadata = re.sub(r"(CHARGE.*\n)", r"\1" + str(comment), metadata)

        # Removing useless "energyN" words
        completed = re.sub("energy0\n", "", completed)
        completed = re.sub("energy1\n", "", completed)
        completed = re.sub("energy2\n", "", completed)
        completed = re.sub("ID=.*\n", "", completed)
        completed = re.sub("\n\n(.*(\n|$)){1,}", "", completed)
        # Removing DONE

        peaks = re.finditer("\n([0-9]{1,20}\.[0-9]{1,10}) ([0-9]{1,4}\.[0-9]{1,20})", completed) # Recover all peaks from spectrum to a peak list.

        dico = {'column1': [], 'column2': []} # Init dictionary of peaks 'mz' 'intensity'

        for matches in peaks: # Saving peaks into dictionary
            dico['column1'].append(float(matches.groups()[0]))
            dico['column2'].append(float(matches.groups()[1]))

        metadata = re.sub(r"(NUM PEAKS: )(.*)", r"NUM PEAKS: " + str(len(dico['column1'])), metadata)  # Write current NUM PEAKS into metadata

        df = pd.DataFrame.from_dict(dico) # Dictionary to pandas dataframe
        df = df.sort_values(by='column1', ascending=0) # Sorting peaks
        dico = df.to_dict('list')

        completed = re.sub(r"(#PMass=[0-9]{1,4}\.[0-9]{1,10})([\s\S]*)", r"\1\n", completed)

        CSV = ""
        for rows in range(len(dico['column1'])):
            CSV = CSV + str(dico['column1'][rows]) + " " + str(dico['column2'][rows]) + "\n"

        completed = re.sub(r"(#PMass=[0-9]{1,4}\.[0-9]{1,10})(\n)", r"\1\n" + CSV + "\n",completed)  # a changer plus tard

        inchikey = None
        if re.search("InChiKey=(.*\n)", completed):
            inchikey = re.search("InChiKey=(.*\n)", completed).group(1)

        if inchikey != None:
            metadata = re.sub(r"(INCHIKEY: )(.*\n)", r"\1" + str(inchikey), metadata)

        completed = metadata + CSV + "\n"

        final_list.append(completed)


        with open("../msp_directory/TEST/InSilico/" + re.sub(".msp", "_InSilico.msp", files), "a",encoding="UTF-8") as file_2:
            file_2.write(final_list[-1])




test_file_path = "../msp_directory/TEST/exp/"


if __name__ == "__main__":  # confirms that the code is under main function
    for files in os.listdir(test_file_path):
        if files.endswith(".msp"):
            with open(os.path.join(test_file_path, files), "r", encoding="UTF-8") as file:
                data = file.read()
                mols_list = data.split("\n\n")

            manager = Manager()
            final_list = manager.list()

            # instantiating process with arguments
            procs = []
            for mols in mols_list:
                # print(name)
                proc = Process(target=cfm_id_treatment, args=(mols,final_list,len(mols_list),files,))
                procs.append(proc)
                proc.start()

            # complete the processes
            for proc in procs:
                proc.join()
