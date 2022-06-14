import os
import sys
import fnmatch
#sys.path.insert(0, '/aloy/home/acomajuncosa/programs/hpc') #CHANGE THIS PATH TO YOUR HPC PATH!
sys.path.insert(0, '/aloy/home/agimeno/Software/hpc') #CHANGE THIS PATH TO YOUR HPC PATH!
from hpc import HPC
from starlife_config import config as cluster_config

scratch_path = "/slgpfs/scratch/irb35/agimeno/MurD/docking/files" 
script_path = "/slgpfs/projects/irb35/agimeno/MurD/docking/files/gen_fps_rDock_center.py"


elements = []
#elements.append(['SagaMurD_Frag373', "/slgpfs/projects/irb35/agimeno/MurD/docking/lib/frags/frags.sdf", 'frags', "/slgpfs/projects/irb35/agimeno/MurD/docking/structures/SagaMurD_Frag373/centroid.sd"])
#elements.append(['SagaMurD_Frag373', "/slgpfs/projects/irb35/agimeno/MurD/docking/lib/subset100/subset100.sdf", 'subset100', "/slgpfs/projects/irb35/agimeno/MurD/docking/structures/SagaMurD_Frag373/centroid.sd"])
for root, dirs, files in os.walk('/slgpfs/projects/irb35/agimeno/MurD/docking/lib', topdown=False):
    for name in files:
        if not fnmatch.fnmatch(name, '2D_' + '*' + '.sdf'): # Do not consider inputs
            if not fnmatch.fnmatch(name, '*' + '_failed.sdf'): # Do not consider failed compounds
                file_path = os.path.join(root, name)
                dir_path = os.path.dirname(file_path)
                dir_name = dir_path.split('/')[-1]
                file_list = [
                    'SagaMurD_Frag373',
                    file_path, 
                    dir_name, 
                    '/slgpfs/projects/irb35/agimeno/MurD/docking/structures/SagaMurD_Frag373/centroid.sd'
                ]
                #print(file_list)
                elements.append(file_list)

ncpus = 1
cluster = HPC(**cluster_config)
njobs = len(elements)

cluster_params = {}
cluster_params['job_name'] = 'murd_docking'
cluster_params["jobdir"] = scratch_path
#cluster_params["memory"] = ncpus
cluster_params['cpu'] = ncpus
cluster_params["wait"] = False
cluster_params["elements"] = elements
cluster_params["num_jobs"] = len(elements)


singularity_image = "/slgpfs/projects/irb35/agimeno/singularity/rDock_image_2/rDock_image_2.simg"
command = "singularity exec {} python {} <TASK_ID> <FILE>".format(
singularity_image,
script_path)

cluster.submitMultiJob(command, **cluster_params)
