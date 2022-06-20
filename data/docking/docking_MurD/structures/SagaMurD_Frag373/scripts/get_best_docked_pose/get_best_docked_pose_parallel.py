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


# # 1. Get the best docked pose for each compound

# Create an sdf file with the best docked pose for each split

# In[ ]:


import fnmatch
import os
import pybel
import subprocess
import re
import time
import math

start_time = time.time()

def sorted_nicely( l ): 
    """ Sort the given iterable in the way that humans expect.""" 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

# Extract best docked pose for each compound
print ('\nGetting the best docked pose obtained for each compound...\n')
results_dir = wd + 'results' + '/'
output_sdf = results_dir + 'results.sd'
output_sdf_best_docked_pose_filename = results_dir + 'best_docked_pose.sdf'
output_sdf_best_docked_pose = pybel.Outputfile('sdf', output_sdf_best_docked_pose_filename, overwrite=True)
if os.path.isfile(output_sdf) == True: # only process the docking results if they exist
    # Generate a dictionary with the best docking score for each compound
    sdf_file = pybel.readfile("sdf",  output_sdf)
    docking_score_dict = {}
    for compound in sdf_file:
        data_dict = compound.data
        compound_title = compound.title
        docking_score = float(data_dict["SCORE.INTER"])
        if compound_title not in docking_score_dict:
            docking_score_dict[compound_title] = docking_score
        else:
            if docking_score <= docking_score_dict[compound_title]:
                docking_score_dict[compound_title] = docking_score
    #print(docking_score_dict)
    # Write the docked pose with the best docking score to the output file
    processed_compounds_list = []
    n=0
    sdf_file = pybel.readfile("sdf",  output_sdf)
    for compound in sdf_file:
        data_dict = compound.data
        compound_title = compound.title
        docking_score = float(data_dict["SCORE.INTER"])
        if compound_title not in processed_compounds_list:
            if math.isclose(docking_score,  docking_score_dict[compound_title]) == True: # compare floats for almost-equality
                output_sdf_best_docked_pose.write(compound)
                n+=1 # counts unique compounds
                print("compound " + str(n) + " of " +str(len(docking_score_dict)) + " written at t = %s seconds" % (time.time() - start_time))
                processed_compounds_list.append(compound_title)

