# conda activate zairachem

import csv
import os
import pandas as pd
from rdkit.Chem import AllChem
from rdkit import Chem
from rdkit.DataStructs import BulkTanimotoSimilarity
import numpy as np
from tqdm import tqdm


RESULTS_FOLDER = "../../results/"
data_folder = "../../data/generated"

# gather all molecules used in the generative model

def read_file(file_name):
    with open(file_name, "r") as f:
        reader = csv.reader(f)
        smiles = []
        for r in reader:
            smiles += [r[0]]
    return smiles

all_used_smiles = []
for file_name in os.listdir(data_folder):
    if file_name.endswith(".csv"):
        all_used_smiles += read_file(os.path.join(data_folder, file_name))
all_used_smiles = list(set(all_used_smiles))
all_used_mols = [Chem.MolFromSmiles(smi) for smi in tqdm(all_used_smiles)]

def calc_fingerprint(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius=3)

all_used_fps = [calc_fingerprint(mol) for mol in tqdm(all_used_mols)]

# other interesting subsets

known_fps = [calc_fingerprint(Chem.MolFromSmiles(smi)) for smi in read_file("../../data/generated/known_hits.csv")]
pocketvec_fps = [calc_fingerprint(Chem.MolFromSmiles(smi)) for smi in read_file("../../data/generated/pocketvec_hits.csv")]

# do the screening

data_folder = os.path.join(data_folder, "reinvent")

for input_file in os.listdir(data_folder):
    if not input_file.endswith(".csv"): continue
    output_file = os.path.join(RESULTS_FOLDER, "tanimoto-{0}".format(input_file))
    if os.path.exists(output_file):
        continue
    with open(os.path.join(data_folder, input_file), "r") as f:
        reader = csv.reader(f)
        smiles = []
        for r in reader:
            smiles += [r[0]]
    mols = [Chem.MolFromSmiles(smi) for smi in smiles]
    R = []
    for mol, smi in tqdm(zip(mols, smiles)):
        if mol is not None:
            ik = Chem.MolToInchiKey(mol)
            fp = calc_fingerprint(mol)
            R += [[ik, smi, np.max(BulkTanimotoSimilarity(fp, all_used_fps)), np.max(BulkTanimotoSimilarity(fp, known_fps)), np.max(BulkTanimotoSimilarity(fp, pocketvec_fps))]]
        else:
            R += [["", "", "", "", ""]]
    df = pd.DataFrame(R, columns=["key", "input", "ts_train", "ts_known", "ts_pocketvec"])
    df.to_csv(output_file, index=False, sep=",")
