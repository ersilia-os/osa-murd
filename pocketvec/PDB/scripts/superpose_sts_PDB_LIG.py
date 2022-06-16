from Bio.PDB import *
from Bio.PDB import PDBList
from Bio.SeqUtils import seq1
import tarfile
import numpy as np
import pickle
import time
import sys
import shutil
import os
import Bio
import tarfile
import gzip
from rdkit.Chem import rdFreeSASA
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from collections import Counter
import pandas as pd
import urllib
from lxml import etree
import requests as r
from Bio import SeqIO
from io import StringIO
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances


task_id = sys.argv[1]  # <TASK_ID>
filename = sys.argv[2]  # <FILE>  

input_pickle = pickle.load(open(filename, 'rb'))

uniprot = input_pickle[task_id][0][0]  # uniprot
structures = input_pickle[task_id][0][1]  # associated structures
ref_structure = input_pickle[task_id][0][2]  # ref structure to make superpositions

# Some functions

# Read trans/rot matrix for alignment (TM-align)
def read_mtx(file):
    f = open(file, "r").readlines()[2:5]
    f = np.array([i.split()[1:] for i in f], dtype='float32')
    t = np.array(f[:,0])
    u = np.array([f[0][1:], f[1][1:], f[2][1:]])
    return t, u

# Get 3D points from ligand centroid
def get_3d_points(file):
    point = PDBParser().get_structure("point", file)[0]
    return [at.get_coord() for at in point.get_atoms()][0]

# Get line to merge pdbs
def get_line_point(file):
    f = open(file, "r").readlines()
    return f[0]






### 1. CREATE DIRECTORIES

path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/structures")


mtx_path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/alignment/MSA/mtx", uniprot)
outpath = os.path.join('/aloy/home/acomajuncosa/MurD/GitHub/alignment/pockets', uniprot)
vis_path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/alignment/MSA/visualization", uniprot)
path_to_aligned_sts = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/alignment/MSA/aligned_structures/", uniprot)
path_to_reports = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/alignment/MSA/reports/", uniprot)

if os.path.exists(mtx_path) is False: os.makedirs(mtx_path)
if os.path.exists(outpath) is False: os.makedirs(outpath)
if os.path.exists(vis_path) is False: os.makedirs(vis_path)
if os.path.exists(path_to_aligned_sts) is False: os.makedirs(path_to_aligned_sts)
if os.path.exists(path_to_reports) is False: os.makedirs(path_to_reports)



### 2. EXTRACT & ORGANIZE THE DATA

# Extract the data (st, lig, ctr)
for structure in structures:
    
    pdb = structure.split("_")[0]
    chain = structure.split("_")[1]
    ligand = structure.split("_")[2]
    count = structure.split("_")[3]


    label = "_".join([pdb, chain, ligand])
    
    
    st = label + "/" + label + ".pdb"
    lig = label + "/" + ligand + "_" + count + ".pdb"
    ctr = label + "/" + ligand + "_" + count + "_centroid.pdb"

    tar = tarfile.open(os.path.join(path, pdb[1:3], label + ".tar.gz"), "r")
    tar.extract(st, path=os.path.join(outpath))
    tar.extract(lig, path=os.path.join(outpath))
    tar.extract(ctr, path=os.path.join(outpath))
    tar.close()
    
# Select one st as reference (st2)
shutil.copyfile(os.path.join(outpath, ref_structure, ref_structure + ".pdb"), os.path.join(vis_path, "cluster_REFERENCE.pdb"))


### 3. MAKE STRUCTURAL SUPERPOSITIONS

# Structural superposition against ref (st2)
for structure in structures:
    
    pdb = structure.split("_")[0]
    chain = structure.split("_")[1]
    ligand = structure.split("_")[2]
    count = structure.split("_")[3]
    
    st1 = "_".join(structure.split("_")[:-1])
    count = structure.split("_")[-1]
    
    # Superpose using TM-align
    command = '/aloy/home/acomajuncosa/programs/TM_align/TMalign ' + os.path.join(outpath, st1, st1 + '.pdb') + ' ' + os.path.join(outpath, ref_structure, ref_structure + '.pdb') + ' -m ' + os.path.join(mtx_path, ref_structure + "_" + st1 + ".txt")
    os.system(command)
    

    t, u = read_mtx(os.path.join(mtx_path, ref_structure + "_" + st1 + ".txt"))


    # Protein
    structure = PDBParser().get_structure("st", os.path.join(outpath, st1, st1 + ".pdb"))[0]
    for atom in structure.get_atoms():
        #atom.transform(rotation_matrix, translation_matrix)
        x0, y0, z0 = atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
        x = t[0] + u[0][0]*x0 + u[0][1]*y0 + u[0][2]*z0
        y = t[1] + u[1][0]*x0 + u[1][1]*y0 + u[1][2]*z0
        z = t[2] + u[2][0]*x0 + u[2][1]*y0 + u[2][2]*z0
        atom.coord = np.array([x, y, z])
    io = PDBIO()
    io.set_structure(structure)
    io.save(os.path.join(outpath, st1, st1 + '_aligned.pdb'))

    label = "_".join([pdb, chain, ligand])
    lig = ligand + "_" + count
    ctr = ligand + "_" + count + "_centroid"

    shutil.copyfile(os.path.join(outpath, st1, st1 + '_aligned.pdb'), os.path.join(path_to_aligned_sts, st1 + '_aligned.pdb'))
    
    
    # Centroid
    structure = PDBParser().get_structure("st", os.path.join(outpath, st1, ctr + ".pdb"))
    for atom in structure.get_atoms():
        # atom.transform(rotation_matrix, translation_matrix)
        x0, y0, z0 = atom.get_coord()[0], atom.get_coord()[1], atom.get_coord()[2]
        x = t[0] + u[0][0]*x0 + u[0][1]*y0 + u[0][2]*z0
        y = t[1] + u[1][0]*x0 + u[1][1]*y0 + u[1][2]*z0
        z = t[2] + u[2][0]*x0 + u[2][1]*y0 + u[2][2]*z0
        atom.coord = np.array([x, y, z])
    io = PDBIO()
    io.set_structure(structure)
    io.save(os.path.join(outpath, st1, ctr + '_aligned.pdb'))
    shutil.copyfile(os.path.join(outpath, st1, ctr + '_aligned.pdb'), os.path.join(vis_path, "CTR_" + st1 + "_" + "_".join(ctr.split("_")[1:]) + ".pdb"))



### 4. READ THE COORDINATES FROM ALL CENTROIDS 

labs, coords = [], []

for file in sorted(os.listdir(vis_path)):
    if "centroid" in file and "REFERENCE" not in file:
        labs.append(file)
        coord = get_3d_points(os.path.join(vis_path, file))
        coords.append([coord[0], coord[1], coord[2]])

coords = np.array(coords)


### 5. CLUSTERING

if len(coords) > 1:

    # Clustering
    clustering = AgglomerativeClustering(distance_threshold=12, n_clusters=None,
     affinity='euclidean', linkage='complete').fit(coords)

    # Pretty print
    for cluster in sorted(set(clustering.labels_)):
        with open(os.path.join(vis_path, "cluster_" + str(cluster) + ".pdb"), "w") as f:
            f.write("HEADER\n")
            for lab, numb in zip(labs, clustering.labels_):
                if numb == cluster:
                    l = get_line_point(os.path.join(vis_path, lab))
                    f.write(l)
            f.write("END")

    # Print report(s)
    with open(os.path.join(path_to_reports, uniprot + ".tsv"), "w") as f:
        f.write("\n".join(["\t".join([i,str(j)]) for i,j in zip(labs, clustering.labels_)]))

    report = {i:j for i,j in zip(labs, clustering.labels_)}
    pickle.dump(report, open(os.path.join(path_to_reports, uniprot + ".pkl"), "wb"))


else:  # In case only one st in this uniprot

    # Pretty print
    with open(os.path.join(vis_path, "cluster_" + str(0) + ".pdb"), "w") as f:
        f.write("HEADER\n")
        for lab, numb in zip(labs, [0]):
            l = get_line_point(os.path.join(vis_path, lab))
            f.write(l)
        f.write("END")

    # Print report(s)
    with open(os.path.join(path_to_reports, uniprot + ".tsv"), "w") as f:
        f.write("\n".join(["\t".join([i,str(j)]) for i,j in zip(labs, [0])]))

    report = {i:j for i,j in zip(labs,[0])}
    pickle.dump(report, open(os.path.join(path_to_reports, uniprot + ".pkl"), "wb"))





### 6. CHECK FOR CLASHES 

# Get the nonredundant set of pairs ==> (A, B, C, D): (A,B), (A,C), (A,D), (B,C), (B,D), (C,D)
def get_nonredundant_pairs(ctr):
    pairs = []
    for i in ctrs:
        for j in ctrs:
            if i != j and " ".join([i, j]) not in pairs and " ".join([j, i]) not in pairs:
                pairs.append(" ".join([i, j]))
    return pairs


def check_clash(points, locations, minimum=2):
    dists = np.array([min(i) for i in pairwise_distances(points, locations, metric='euclidean')])
    if min(dists) < minimum:
        return True
    else:
        return False


def check_conflict(loc1, loc2, path_to_st):
    
    # Load loc from all atoms in st
    st = PDBParser().get_structure("st", path_to_st)
    atoms_prot = [i.coord for i in st.get_atoms()]
    
    # Calculate points to check
    thr = 1  # point every 1 A
    incr = np.linalg.norm(loc1-loc2)
    n_spaces = int(incr/thr) + 2
    points_to_check = np.linspace(loc1, loc2, n_spaces)
    
    # Check clashes with points
    conflict = check_clash(points_to_check, atoms_prot)
    
    return conflict


clusters = pickle.load(open(os.path.join(path_to_reports, uniprot + ".pkl"), "rb"))
numb_clusters = sorted(set([clusters[i] for i in clusters]))
report_conflicts = {}

for numb_cluster in numb_clusters:
    ctrs, loc = [], {} 
    
    # Get how many centroids are clustered in numb_clustered
    for ctr in clusters:
        if clusters[ctr] == numb_cluster:
            ctrs.append(ctr)
            
    # Get all pairs of distances
    pairs = get_nonredundant_pairs(ctr)
    for ctr in ctrs:
        loc[ctr] = get_3d_points(os.path.join(vis_path, ctr))
        
    # Check conflicts for all pairs of centroids
    c = False
    for pair in pairs:
        ctr1, ctr2 = loc[pair.split()[0]], loc[pair.split()[1]]
        conflict = check_conflict(ctr1, ctr2, os.path.join(vis_path, "cluster_REFERENCE.pdb"))
        if conflict is True:
            c = True
            break
    report_conflicts[numb_cluster] = [c, len(ctrs)]

pickle.dump(report_conflicts, open(os.path.join(path_to_reports, uniprot + "_conflicts.pkl"), "wb"))