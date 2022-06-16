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
from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt
from collections import Counter
import pandas as pd
import urllib
from lxml import etree
import requests as r
from Bio import SeqIO
from io import StringIO
import signal
import tarfile
from subprocess import Popen, PIPE


task_id = sys.argv[1]  # <TASK_ID>
filename = sys.argv[2]  # <FILE>  

input_pickle = pickle.load(open(filename, 'rb'))

infile = input_pickle[task_id][0][0]  # path_in
prepared_prot = input_pickle[task_id][0][1]  # path_out
logfile = input_pickle[task_id][0][2]  # path_log


# 5. Prepare st MOE
def prepare(infile, prepared_prot, logfile, all_outputs):
    # if os.path.isfile(prepared_prot):
    #     return True, all_outputs
    # else:
    moefunct = "/aloy/home/acomajuncosa/PocketVec_v2/pocketvec/ProSPECCTs/moefunctions.svl"
    # PROTEINS ARE PREPARED USING MOE
    process=Popen(["/aloy/home/acomajuncosa/programs/MOE/moe2020/bin/moebatch -load " + moefunct + " -exec \"proteinprep['" + infile + "','" + prepared_prot + "','" + logfile + "']\""],stdout=PIPE,stderr=PIPE,shell=True)
    stdout, stderr = process.communicate()
    sys.stderr.write(str(stdout) + "\n\n")
    sys.stderr.write(str(stderr) + "\n\n")
    all_outputs.append(str(stderr))
    if os.path.isfile(prepared_prot):
        return True, all_outputs
    else:
        return False, all_outputs

def handler(signum, frame):
    raise Exception("Protein preparation exceeds 15min. Cancelling job...")

#### STRUCTURE PREPARATION ####

prepared = False
timeout = False
all_outputs = []

while prepared is False and timeout is False:

    # Register the signal function handler
    signal.signal(signal.SIGALRM, handler)

    # Define a timeout for your function ==> 1h
    signal.alarm(1000)  ## CHANGE THIS IF A PDB IS TAKING TOO MUCH

    try:
        result, all_outputs = prepare(infile, prepared_prot, logfile, all_outputs)
    except Exception as exc:
        result, timeout = False, True
        sys.stderr.write(str(exc))
        sys.stderr.write("\n\n\n")
        sys.stderr.write("TIMEOUT!\n\n")
        sys.stderr.flush()

    if result is True:
        prepared = True
    elif timeout is False:
        sys.stderr.write("Preparation did not take place as expected. Issues with licenses are the most probable cause. Trying again...\n\n")
        sys.stderr.flush()
        time.sleep(30)