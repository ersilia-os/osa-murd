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


# # 1. Calculate accessibility

# Create an sdf file with the best docked pose for each split

# In[ ]:


import os
import pybel
import time
from rdkit import Chem
from rdkit.Chem import rdFreeSASA
import subprocess
import numpy as np

start_time = time.time()

def SASA(prot, lig): 

    # Protonation gives too many issues. Avoid it
    
    #compute ligand SASA
    #lig_h = Chem.rdmolops.AddHs(lig, addCoords=True, explicitOnly=True)
    lig_h = lig
    # Get Van der Waals radii (angstrom)
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in lig_h.GetAtoms()]
    # Compute solvent accessible surface area
    lig_sasa = rdFreeSASA.CalcSASA(lig_h, radii)

    # Join protein & ligand
    comp = Chem.CombineMols(prot, lig)
    comp_h = comp
    #comp_h = Chem.AddHs(comp, addCoords=True)
    # Get Van der Waals radii (angstrom)
    ptable = Chem.GetPeriodicTable()
    radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in comp_h.GetAtoms()]
    # Compute solvent accessible surface area
    comp_sasa = rdFreeSASA.CalcSASA(comp_h, radii)
    comp_lig = Chem.GetMolFrags(comp_h, asMols=True,  sanitizeFrags=True)
    comp_lig = [i for i in comp_lig if lig_h.GetNumHeavyAtoms() == i.GetNumHeavyAtoms()][0]
    
    lig_sasa_free = 0
    for a in lig_h.GetAtoms():
        lig_sasa_free += float(a.GetProp("SASA"))

    lig_sasa_bound = 0
    for a in comp_lig.GetAtoms():
        lig_sasa_bound += float(a.GetProp("SASA"))
        
    return round(lig_sasa_free, 3), round(lig_sasa_bound, 3)

# Get SASA
#prot = Chem.MolFromPDBFile(os.path.join(path, pdb + "_" + chain_id + "_" + int_lig + ".pdb"))

# # Convert protein from mol2 to pdb
# prot_mol2_file = parent_dir + 'SagaMurD_Frag373_MOEprep.mol2'
# prot_pdb_file = parent_dir + 'SagaMurD_Frag373_MOEprep.pdb'
# obabel_command = 'obabel -imol2 ' + prot_mol2_file + ' -opdb -O ' + prot_pdb_file
# subprocess.call(obabel_command.split())
# prot = Chem.MolFromPDBFile(prot_pdb_file)

# Read prot
prot_pdb_file = parent_dir + 'SagaMurD_Frag373_MOEprep.pdb'
prot = Chem.MolFromPDBFile(prot_pdb_file)

results_dir = wd + 'results' + '/'
# Create dataframe to store the results
df_accessibility = pd.DataFrame(columns=['title','score','accessibility'])

docked_poses_file = results_dir + 'best_docked_pose.sdf'
sdf_file = pybel.readfile("sdf",  docked_poses_file)
for compound in sdf_file:
    row_dict = {}
    data_dict = compound.data
    compound_title = compound.title
    docking_score = float(data_dict["SCORE.INTER"])
    # Convert compound to rdkit mol
    compound.write('sdf', results_dir + 'mol.sdf', overwrite=True)
    suppl = Chem.SDMolSupplier(results_dir + 'mol.sdf')
    mol_list = [x for x in suppl]
    lig = mol_list[0]
    # Convert rdkit mol to pdb
    #lig = Chem.rdmolfiles.MolToPDBBlock(lig)
    
    # Calculate accessibility
    print('prot', prot)
    print('type(prot)', type(prot))
    print('lig', lig)
    print('type(lig)', type(lig))
    
    lig_sasa_free, lig_sasa_bound = SASA(prot, lig)
    acc = round(lig_sasa_bound / lig_sasa_free, 3)

#     try:
#         lig_sasa_free, lig_sasa_bound = SASA(prot, lig)
#         acc = round(lig_sasa_bound / lig_sasa_free, 3)
#     except: # Boost.Python.ArgumentError
#         acc = np.nan
        
    #print('compound.title', 'docking_score', 'acc')
    #print(compound.title, docking_score, acc)
    
    # Append rows to dataframe
    row_dict['title'] = compound.title
    row_dict['score'] = docking_score
    row_dict['accessibility'] = acc
    df_accessibility = df_accessibility.append(row_dict, ignore_index=True)

# Export dataframe as .tsv
df_accessibility_file = results_dir + 'accessibilities' + '.tsv'
df_accessibility.to_csv(df_accessibility_file, sep='\t', encoding='utf-8', index=False)
print(df_accessibility)
    

