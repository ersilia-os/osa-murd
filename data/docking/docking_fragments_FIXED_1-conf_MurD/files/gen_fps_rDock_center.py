import numpy as np
import subprocess
import sklearn
import pickle
import time
import sys
import os
import shutil
import tarfile
import pybel
from collections import Counter
import scipy.stats as ss


task_id = sys.argv[1]
filename = sys.argv[2]

input_pickle = pickle.load(open(filename, 'rb'))

st = input_pickle[task_id][0][0]  # st
path_to_lib = input_pickle[task_id][0][1]  # library to dock. should be in sdf format
lib = input_pickle[task_id][0][2]  # library to dock. name (e.g. LLM, Frags, test, etc)
path_to_ctr = input_pickle[task_id][0][3]  # path to centroid. should be in sd format


radius = "12.0"  # Optimal value. If needed, try increasing it
nruns = 25
path_to_data = "/slgpfs/projects/irb35/agimeno/MurD/docking/structures"
path_to_st = os.path.join(path_to_data, st)
path_to_out = os.path.join(path_to_data, st, 'rDock_results_' + lib, "results")

sys.stderr.write("\n\n\n" + str(st) + "\n\n\n")
sys.stderr.flush()

try:
    if os.path.exists(path_to_out) is False: os.makedirs(path_to_out)
except:
    pass


os.chdir(path_to_out)
os.environ["RBT_ROOT"] = "/opt/rdock"
# os.environ["LD_LIBRARY_PATH"] = "$LD_LIBRARY_PATH:/opt/rdock/lib"
# os.environ["PATH"] = "$PATH:/opt/rdock/bin"

sys.stderr.write("\n\n" + str(os.getenv('RBT_ROOT')) + "\n\n")
sys.stderr.write("\n\n" + str(os.getenv('RBT_HOME')) + "\n\n")
sys.stderr.write("\n\n" + str(os.getenv('LD_LIBRARY_PATH')) + "\n\n")
sys.stderr.write("\n\n" + str(os.getenv('PATH')) + "\n\n")
sys.stderr.flush()

################################
#---The real script starts here
################################

# 1. Generate structure parameter file

text = """RBT_PARAMETER_FILE_V1.00
TITLE Change this if you want

RECEPTOR_FILE """ + os.path.join(path_to_st, st + "_MOEprep.mol2") + """
RECEPTOR_FLEX 3.0

#################################################################
## CAVITY DEFINITION: REFERENCE LIGAND METHOD
#################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL """ + path_to_ctr + """
    RADIUS """ + radius + """
    SMALL_SPHERE 1.0
    MIN_VOLUME 100
    MAX_CAVITIES 1
END_SECTION


SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION

#################################################
# WARNING!!! WE ARE DOING RIGID LIGAND DOCKING!!
# CHANGE DIHEDRAL MODE "FIXED" TO "FREE" IN 
# ORDER TO PERFORM FLEXIBLE DOCKING
#################################################
SECTION LIGAND
    TRANS_MODE FREE
    ROT_MODE FREE
    DIHEDRAL_MODE """ + "FIXED" + """
END_SECTION"""


with open(path_to_out + "/receptor_parameters.prm", "w") as f:
    f.write(text)


# 2. Define cavity  ---->  generated files: .as & .grd
command = 'rbcavity -was -d -r ./receptor_parameters.prm > ./cavity_log.log'
os.system(command)

os.rename(os.path.join(path_to_out, "_cav1.grd"), os.path.join(path_to_out, "cavity_" + st + ".grd"))


# 3. Run rDock using std parameters
shutil.copy("/slgpfs/projects/irb35/acomajuncosa/PocketVec_v2/data/dock.prm", path_to_out + "/dock.prm")


command = 'rbdock -i ' + path_to_lib + ' -o ./results -r ./receptor_parameters.prm -p ./dock.prm -n ' + str(nruns) + ' -s 42 > ./rdock_log.log'
os.system(command)








with open(os.path.join(path_to_st, 'rDock_results_' + lib, 'scores'), "w") as out_score:

    # Look for errors
    errors = False
    with open(os.path.join(path_to_out, "rdock_log.log"), "r") as f:
        for l in f:
            if "error" in l.lower() or ("warning" in l.lower() and "undefined tripos type" not in l.lower()):
                errors = True

    # If no errors, select affinity
    if errors is False:
        file = pybel.readfile("sd", os.path.join(path_to_out, "results.sd"))
        aff = 1000000
        affinities = {}
        for molecule in file:
            if molecule.title not in affinities:
                affinities[molecule.title] = aff
            value = float(molecule.data['SCORE.INTER'])
            if value < affinities[molecule.title]:
                affinities[molecule.title] = value
        
        for molecule in sorted(affinities):
            out_score.write(molecule + "\t" + str(affinities[molecule]) + "\n")


# Tar results
os.chdir(os.path.join(path_to_st, 'rDock_results_' + lib))
command = "tar -czf results.tar.gz results"
os.system(command)

# Remove folders
shutil.rmtree("results")
