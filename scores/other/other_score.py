# conda activate zairachem

import os
import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.QED import qed
from rdkit.Chem import Descriptors

RESULTS_FOLDER = "../../results/"
data_folder = "../../data/generated/reinvent"

for input_file in os.listdir(data_folder):
    if not input_file.endswith(".csv"): continue
    output_file = os.path.join(RESULTS_FOLDER, "other-{0}".format(input_file))
    if os.path.exists(output_file):
        continue
    with open(os.path.join(data_folder, input_file), "r") as f:
        reader = csv.reader(f)
        smiles = []
        for r in reader:
            smiles += [r[0]]
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    R = []
    for mol, smi in zip(mols, smiles):
        if mol is not None:
            ik = Chem.MolToInchiKey(mol)
            R += [[ik, smi, Descriptors.MolWt(mol), MolLogP(mol), qed(mol)]]
        else:
            R += [["", "", "", "", ""]]
    df = pd.DataFrame(R, columns=["key", "input", "mw", "clogp", "qed"])
    df.to_csv(output_file, index=False, sep=",")