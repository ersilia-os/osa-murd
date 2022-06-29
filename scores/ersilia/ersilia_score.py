# conda activate ersilia

import os
from ersilia import ErsiliaModel

from rdkit import Chem
import csv

RESULTS_FOLDER = "../../results/"

model_ids = ["eos2r5a", "eos7pw8"]

data_folder = "../../data/generated/intermediate"

for input_file in os.listdir(data_folder):
    if not input_file.endswith(".csv"): continue
    for model_id in model_ids:
        output_file = os.path.abspath(os.path.join(RESULTS_FOLDER, "ersilia-{0}-{1}".format(model_id, os.path.split(input_file)[-1])))
        input_file = os.path.abspath(os.path.join(data_folder, input_file))
        if os.path.exists(output_file):
            continue
        print(input_file, output_file)
        with open(input_file, "r") as f:
            smiles = []
            reader = csv.reader(f)
            for r in reader:
                mol = Chem.MolFromSmiles(r[0])
                if mol is not None:
                    smiles += [r[0]]
                else:
                    smiles += ["CCCCOCCCC"]
        with open("tmp.csv", "w") as f:
            writer = csv.writer(f)
            for smi in smiles:
                writer.writerow([smi])
        em = ErsiliaModel(model=model_id)
        em.serve()
        em.api(input="tmp.csv", output=output_file)
        em.close()


os.remove("tmp.csv")
