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
labs = np.array(labs)


### 5. CLUSTERING

# Some functions for clustering

def check_cond2(cluster_points, max_centroid):
    # Get new centroid
    ctr = calculate_centroid(cluster_points)
    dist = max(pairwise_distances(ctr.reshape(1, -1), cluster_points)[0])
    return dist < max_centroid

def get_min_dist_between_clusters(points1, points2):
    return min([j for i in pairwise_distances(points1, points2) for j in i])

def calculate_centroid(arr):
    arr = np.array(arr)
    length = arr.shape[0]
    sum_x = np.sum(arr[:, 0])
    sum_y = np.sum(arr[:, 1])
    sum_z = np.sum(arr[:, 2])
    return np.array([sum_x/length, sum_y/length, sum_z/length])


def clustering(coords, min_newpoint, max_centroid):

    # Get info
    coords = {i: j for i, j in zip(range(len(coords)), coords)}
    # clusters = {i: j.reshape(1, -1) for i, j in zip(range(len(coords)), coords)}
    real_clusters = {}
    clust_gen = True
    clust_count = -1

    while clust_gen and len([i for i in coords if i not in real_clusters]) > 1:

        # Get minimum distance between points
        labels = [i for i in coords if i not in real_clusters]
        distances = pairwise_distances([coords[i] for i in coords if i not in real_clusters])
        mindist = {}
        for c1, row in zip(labels, distances):
            for c2, dist in zip(labels, row):
                c1, c2 = str(c1), str(c2)
                if c1 + "_" + c2 not in mindist and c2 + "_" + c1 not in mindist and c1 != c2:
                    mindist[c1 + "_" + c2] = dist

        # Get the first centroid
        minimum = sorted(mindist, key=lambda x: mindist[x])[0]
        # print(mindist[minimum], clust_count+1)
        # print(clust_gen, cont)

        if mindist[minimum] < min_newpoint:
            #print("in!")
            clust_count += 1
            real_clusters[int(minimum.split("_")[0])] = clust_count
            real_clusters[int(minimum.split("_")[1])] = clust_count
            if len([i for i in coords if i not in real_clusters]) > 0:
                clust_gen = True
                cont = True
            else:
                clust_gen = False
                cont = False

        else:
            clust_gen = False
            cont = False

        while cont == True and len([i for i in coords if i not in real_clusters]) > 0:

            cluster_points = np.array([coords[i] for i in real_clusters if real_clusters[i] == clust_count])
            cluster_labels = np.array([i for i in real_clusters if real_clusters[i] == clust_count])
            centroid = calculate_centroid(cluster_points)

            ### Increase cluster
            # Get the min external point to cluster points
            increasing = {}
            labels = np.array([i for i in coords if i not in real_clusters])
            for cluster_label, cluster_point in zip(cluster_labels, cluster_points):
                dist = pairwise_distances(cluster_point.reshape(1, -1), [coords[i] for i in coords if i not in real_clusters])[0]
                ind = dist.argmin()
                increasing[cluster_label] = [dist[ind], labels[ind]]


            increasing = np.array([increasing[i] for i in increasing])
            dist, lab = increasing[:,0], increasing[:,1]
            ind = dist.argmin()
            dist, lab = dist[ind], int(lab[ind])

            # Check condition 1 // cond1
            if dist < min_newpoint:
                cond1 = True
            else:
                cond1 = False


            # Check condition 2 // cond2
            if cond1 == True:
                cluster_points_check = np.array([list(coords[i]) for i in real_clusters if real_clusters[i] == clust_count] + [list(coords[lab])])
                cond2 = check_cond2(cluster_points_check, max_centroid)

            if cond1 == True and cond2 == True:
                real_clusters[lab] = clust_count
            else:
                cont = False

    return real_clusters, clust_count


def assign_outliers(clust_count, real_clusters, coords):        

    clust_count += 1
    coords = {i: j for i, j in zip(range(len(coords)), coords)}

    for i in sorted(coords):
        if i not in real_clusters:
            real_clusters[i] = clust_count
            clust_count += 1
            
    return real_clusters


def pretty_clusters(real_clusters, vis_path, path_to_reports):
    
    clst = {}

    for numb_cluster in set(sorted([real_clusters[i] for i in real_clusters])):
        for st in real_clusters:
            if real_clusters[st] == numb_cluster:
                clst[labs[st]] = numb_cluster
                
    # Pretty print
    for cluster in sorted(set([real_clusters[i] for i in real_clusters])):
        with open(os.path.join(vis_path, "cluster_" + str(cluster) + ".pdb"), "w") as f:
            f.write("HEADER\n")
            for lab in clst:
                if clst[lab] == cluster:
                    l = get_line_point(os.path.join(vis_path, lab))
                    f.write(l)
            f.write("END")

    # Print report(s)
    with open(os.path.join(path_to_reports, uniprot + ".tsv"), "w") as f:
        f.write("\n".join(["\t".join([i,str(clst[i])]) for i in clst]))

    report = {i:clst[i] for i in clst}
    pickle.dump(report, open(os.path.join(path_to_reports, uniprot + ".pkl"), "wb"))


# Actual clustering

min_newpoint = 5
max_centroid = 18

if len(coords) > 1:

    real_clusters, clust_count = clustering(coords, min_newpoint, max_centroid)
    real_clusters = assign_outliers(clust_count, real_clusters, coords)
    pretty_clusters(real_clusters, vis_path, path_to_reports)

else:

    real_clusters = {}
    real_clusters[0] = 0
    pretty_clusters(real_clusters, vis_path, path_to_reports)