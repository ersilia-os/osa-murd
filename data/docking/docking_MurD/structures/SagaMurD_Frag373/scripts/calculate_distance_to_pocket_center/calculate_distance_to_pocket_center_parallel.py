#!/usr/bin/env python
# coding: utf-8

# # 0. Definitions

# In[ ]:


import os
import pandas as pd
import numpy as np
import random
import pickle

# get the notebook's root path
try: 
    ipynb_path
except NameError: 
    ipynb_path = os.getcwd()
if ipynb_path.startswith('/slgpfs/'):
    # change scratch path to projects path
    ipynb_path = ipynb_path.replace("/scratch/", "/projects/")
    #print(ipynb_path)
elif ipynb_path.startswith('/aloy/'):
    # change scratch path to home path
    ipynb_path = ipynb_path.replace("/scratch/", "/home/")
    #print(ipynb_path)
parent_dir = ipynb_path + '/' + '../' + '../'
#print('wd:', wd)


# Get dirname

# In[ ]:


import sys
import time
import json

task_id = sys.argv[1]  # <TASK_ID>
filename = sys.argv[2]  # <FILE>
input_pickle = pickle.load(open(filename, 'rb'))
element = input_pickle[task_id][0] #This "[0]" is important as it by default uses a list of lists

# Get dirname
wd = parent_dir + element + '/'


# # 1. Calculate distance to pocket center

# Calculate distance to pocket center

# In[ ]:


import os
import pybel
import time
from rdkit import Chem
from rdkit.Chem import rdFreeSASA
import subprocess
import numpy as np

start_time = time.time()

# # Read prot
# prot_pdb_file = parent_dir + 'SagaMurD_Frag373_MOEprep.pdb'
# prot = Chem.MolFromPDBFile(prot_pdb_file)

# Obtain the coordinates of the center of the pocket (use the center of the grid)
pocket_coords = np.array([32.729, -33.675, -1.107])

results_dir = wd + 'results' + '/'
# Create dataframe to store the results
df_pocket_edist = pd.DataFrame(columns=['title','score','pocket_edist'])

docked_poses_file = results_dir + 'best_docked_pose.sdf'
sdf_file = pybel.readfile("sdf",  docked_poses_file)
for compound in sdf_file:
    row_dict = {}
    data_dict = compound.data
    compound_title = compound.title
    docking_score = float(data_dict["SCORE.INTER"])
#     # Convert compound to rdkit mol
#     compound.write('sdf', results_dir + 'mol.sdf', overwrite=True)
#     suppl = Chem.SDMolSupplier(results_dir + 'mol.sdf')
#     mol_list = [x for x in suppl]
#     lig = mol_list[0]
    # Convert rdkit mol to pdb
    #lig = Chem.rdmolfiles.MolToPDBBlock(lig)
    atom_edists_dict = {}
    # Calculate euclidean distance to the center of the pocket for each atom
    for atom in compound:
        atom_idx = str(atom.idx)
        atom_coords = atom.coords
        #print('atom_coords', atom_coords)
        #print('type(atom_coords)', type(atom_coords))
        atom_edist = np.linalg.norm(pocket_coords-atom_coords)
        atom_edists_dict[atom_idx] = atom_edist
    # Get the minimum distance
    key_min = min(atom_edists_dict.keys(), key=(lambda k: atom_edists_dict[k]))
    min_edist = atom_edists_dict[key_min]
    pocket_edist = min_edist
    #print('compound.title', 'docking_score', 'pocket_edist')
    #print(compound.title, docking_score, pocket_edist)
    
    # Append rows to dataframe
    row_dict['title'] = compound.title
    row_dict['score'] = docking_score
    row_dict['pocket_edist'] = pocket_edist
    df_pocket_edist = df_pocket_edist.append(row_dict, ignore_index=True)

# Export dataframe as .tsv
df_pocket_edist_file = results_dir + 'pocket_edist' + '.tsv'
df_pocket_edist.to_csv(df_pocket_edist_file, sep='\t', encoding='utf-8', index=False)
print(df_pocket_edist)
    

