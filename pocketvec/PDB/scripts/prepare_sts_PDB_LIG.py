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

#sys.stderr.write(Bio.__version__)
#sys.stderr.flush()


#time.sleep(float(random.randint(0, 100))/10)

task_id = sys.argv[1]  # <TASK_ID>
filename = sys.argv[2]  # <FILE>  

input_pickle = pickle.load(open(filename, 'rb'))

infos = input_pickle[task_id][0]  # Information

path_to_summary = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/summary", "summary_" + str(task_id) + ".tsv")


def get_centroid(infile, outfile):
    
    parser = PDBParser()
    structure = parser.get_structure("lig", infile)

    res = [i for i in structure.get_residues()][0]
    ligatoms = [at.coord for at in res.get_atoms()]
    x = np.mean(np.array(ligatoms)[:,0])
    y = np.mean(np.array(ligatoms)[:,1])
    z = np.mean(np.array(ligatoms)[:,2])
    center = np.array([x, y, z], dtype=np.float32)

    x, y, z = str(round(x, 3)), str(round(y, 3)), str(round(z, 3))
    ctr = " "*(8-len(x)) + x + " "*(8-len(y)) + y + " "*(8-len(z)) + z
    text = """HEADER\nHETATM    1   C  CTR A   1    """ + ctr + """  1.00  1.00           C\nEND"""

    with open(outfile, "w") as f:
        f.write(text)


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


pdbcode_to_inchi = pd.read_csv("/aloy/home/acomajuncosa/MurD/GitHub/osa-murd/data/PDB/mapping/Components-inchi.ich.txt", sep="\t", header=None, names=['inchi', 'PDB', 'name'], usecols=[0,1])  
# From df to dict's  // PDB-LIGs
d = {}
for i,j in zip(pdbcode_to_inchi['inchi'], pdbcode_to_inchi['PDB']):
    if str(j) == 'nan': j = "NA"
    d[str(j)] = i
pdbcode_to_inchi = d; del d




with open(path_to_summary, "w") as outfile:


    for info in infos:


        try:


            pdb = info[0]
            int_lig = info[1]
            lig = info[2]
            het = info[3]

            path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/structures", pdb[1:3], pdb + "_" + int_lig)
            if os.path.exists(path) is False: os.makedirs(path)
            os.chdir(os.path.join(path))

            ### 1. Create directory & download pdb structure
            with gzip.open('/aloy/home/acomajuncosa/programs/localpdb/mirror/pdb/' + pdb[1:3].lower() + '/pdb' + pdb.lower() + '.ent.gz', 'rb') as f_in:
                with open(os.path.join(path, pdb + ".pdb"), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)


            ligands_to_remove = ['W']
            ligands_to_remove.append("H_" + " "*(3-len(int_lig)) + int_lig)
            ligands_to_remove.append("H_" + int_lig + " "*(3-len(int_lig)))

            if type(lig) == list:
                for l in lig:
                    ligands_to_remove.append("H_" + " "*(3-len(l)) + l)
                    ligands_to_remove.append("H_" + l + " "*(3-len(l)))
                ligands_to_remove = set(ligands_to_remove)


            # 2. Remove water molecules and ligands // but not non std residues!
            parser = PDBParser()
            path_in = os.path.join(path, pdb + ".pdb")
            path_out = path
            structure = parser.get_structure("st", path_in)[0] # Take only the first model. Trivial for X-ray, only 1st for NMR

            # Take only the first model
            io = PDBIO()
            io.set_structure(structure)
            io.save(os.path.join(path, pdb + "_model.pdb"))

            parser = PDBParser()
            path_in = os.path.join(path, pdb + "_model.pdb")
            path_out = path
            structure = parser.get_structure("st", path_in)

            class remove_ligs(Select):
                def accept_residue(self, residue):
                    if residue.get_id()[0] not in ligands_to_remove:
                        return 1
                    else:
                        return 0

            io = PDBIO()
            io.set_structure(structure)
            io.save(os.path.join(path, pdb + "_st.pdb"), remove_ligs())

            # 3. Remove hydrogens
            command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + os.path.join(path_out , pdb + '_st.pdb') + ' -o ' + os.path.join(path_out , pdb + '_hydrogens.pdb') + ' --force_save --non_interactive rem_hydrogen --remove Yes'
            o = os.popen(command).read()
            sys.stderr.write(o + "\n\n")
            sys.stderr.flush()

            # 4. Select occupancies
            command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + os.path.join(path_out , pdb + '_hydrogens.pdb') + ' -o ' + os.path.join(path_out , pdb + '_altloc.pdb') + ' --force_save --non_interactive altloc --select occupancy'
            o = os.popen(command).read()
            sys.stderr.write(o + "\n\n")
            sys.stderr.flush()

            shutil.copyfile(os.path.join(path_out , pdb + '_altloc.pdb'), os.path.join(path_out , pdb + '_' + int_lig + '.pdb'))



            parser = PDBParser()
            structure = parser.get_structure("st", os.path.join(path, pdb + "_model.pdb"))
            interesting_ligands = [i for i in structure.get_residues() if i.get_resname() == int_lig]


            for c, interesting_ligand in enumerate(interesting_ligands):

                try:
                
                    # 5. Select ligand and save it separately
                    class LigSelect(Select):
                        def accept_residue(self, residue):
                            if residue == interesting_ligand:
                                return 1
                            else:
                                return 0
                            
                    io = PDBIO()
                    io.set_structure(structure)
                    io.save(os.path.join(path, int_lig + "_" + str(c) + ".pdb"), LigSelect())
                    #command = 'obabel ' + os.path.join(path, int_lig + "_" + str(c) + ".pdb") + " -O " + os.path.join(path, int_lig + "_" + str(c) + ".sd")
                    #os.system(command)
                    
                    # Get centroid
                    get_centroid(os.path.join(path, int_lig + "_" + str(c) + ".pdb"), os.path.join(path, int_lig + "_" + str(c) + "_centroid.pdb"))
                    command = 'obabel ' + os.path.join(path, int_lig + "_" + str(c) + "_centroid.pdb") + " -O " + os.path.join(path, int_lig + "_" + str(c) + "_centroid.sd")
                    os.system(command)


                    # Get SASA
                    prot = Chem.MolFromPDBFile(os.path.join(path, pdb + "_" + int_lig + ".pdb"))
                    lig = Chem.MolFromPDBFile(os.path.join(path, int_lig + "_" + str(c) + ".pdb"))

                    if lig.GetNumHeavyAtoms() == Chem.MolFromInchi(pdbcode_to_inchi[int_lig]).GetNumHeavyAtoms():
                        
                        lig_sasa_free, lig_sasa_bound = SASA(prot, lig)

                    else:
                        
                        lig_sasa_free, lig_sasa_bound = np.nan, np.nan
                        
                    acc = round(lig_sasa_bound / lig_sasa_free, 3)

                    # Write results
                    outfile.write("\t".join([pdb, int_lig, str(c), interesting_ligand.get_parent().id, str(lig_sasa_free), str(lig_sasa_bound), str(acc)]) + "\n")

                except:

                    # Write results
                    outfile.write("\t".join([pdb, int_lig, str(c), interesting_ligand.get_parent().id, 'failed2', 'failed2', 'failed2']) + "\n")

        except:

            # Failed in the starting steps..
            outfile.write("\t".join([pdb, int_lig, "-", "-", 'failed1', 'failed1', 'failed1']) + "\n")



        path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/structures", pdb[1:3])
        os.chdir(path)   
        tar = tarfile.open(pdb + "_" + int_lig + ".tar.gz", "w:gz")
        tar.add(pdb + "_" + int_lig)
        shutil.rmtree(pdb + "_" + int_lig)
        tar.close()
