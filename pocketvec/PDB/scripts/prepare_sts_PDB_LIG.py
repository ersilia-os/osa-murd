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

#sys.stderr.write(Bio.__version__)
#sys.stderr.flush()


#time.sleep(float(random.randint(0, 100))/10)

task_id = sys.argv[1]  # <TASK_ID>
filename = sys.argv[2]  # <FILE>  

input_pickle = pickle.load(open(filename, 'rb'))

infos = input_pickle[task_id][0]  # Information

path_to_summary = os.path.join("/aloy/home/acomajuncosa/PocketVec_v2/kinase/PDB/LIG/preprocess/summary", "summary_" + str(task_id) + ".tsv")
dict_coverages = pickle.load(open("/aloy/home/acomajuncosa/PocketVec_v2/kinase/PDB/dict_coverages.pkl", "rb"))


with open(path_to_summary, "w") as outfile:


    for info in infos:

        sys.stderr.write(str(info) + "\n\n\n\n")
        sys.stderr.flush()


        pdb = info['PDB']
        uniprot = info['Uniprot']
        pfam = info['Pfam']
        intlig = info['Int Lig']
        chain = info['Chain']
        st_res = int(info['St_res'])
        end_res = int(info['End_res'])

        domain = '_'.join([uniprot, pfam, str(st_res), str(end_res)])

        path = os.path.join("/aloy/home/acomajuncosa/PocketVec_v2/kinase/PDB/LIG/preprocess/data/", domain)

        try:
            if os.path.exists(path) is False: os.makedirs(path)
        except:
            pass


        #####################################
        #-----PARSING PDB & ACCESSIBILITY----
        #####################################


        label = domain + "_" + pdb + "_" + chain + "_" + intlig
        got_pdb = False

        if True:#try:

            ## 1. Create directory & download pdb structure
            if os.path.exists(os.path.join(path, label)) is False: os.makedirs(os.path.join(path, label))
            os.chdir(os.path.join(path, label))

            with gzip.open('/aloy/home/acomajuncosa/programs/localpdb/mirror/pdb/' + pdb[1:3].lower() + '/pdb' + pdb.lower() + '.ent.gz', 'rb') as f_in:
                with open(os.path.join(path, label, pdb + ".pdb"), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            time.sleep(1)
            got_pdb = True


            ## 2.1. Remove waters
            command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + pdb + '.pdb -o ' + pdb + '_water.pdb --force_save --non_interactive water --remove Yes'
            o = os.popen(command).read()
            sys.stderr.write(o + "\n\n\n")
            sys.stderr.flush()


            ## 2.2. Remove hydrogens
            command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + pdb + '_water.pdb -o ' + pdb + '_hydrogens.pdb --force_save --non_interactive rem_hydrogen --remove Yes'
            o = os.popen(command).read()
            sys.stderr.write(o + "\n\n")
            sys.stderr.flush()


            ## 2.3. Select occupancies
            command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + pdb + '_hydrogens.pdb' + ' -o ' + pdb  + '_altloc.pdb --force_save --non_interactive altloc --select occupancy'
            o = os.popen(command).read()
            sys.stderr.write(o + "\n\n")
            sys.stderr.flush()


            ## 3. Select domain
            class DomSelect(Select):
                def accept_residue(self, residue):
                    if residue.get_id()[1] in residues and residue.get_parent().id == chain:
                        return 1
                    else:
                        return 0

            residues = np.array([dict_coverages["_".join([pdb.lower(), uniprot, chain, pfam, str(st_res), str(end_res)])][i] for i in range(st_res, end_res+1) if dict_coverages["_".join([pdb.lower(), uniprot, chain, pfam, str(st_res), str(end_res)])][i][0] is not None])[:,0]
            residues = set([int(i) for i in residues])

            parser = PDBParser()
            structure = parser.get_structure("st", os.path.join(path, label, pdb + "_altloc.pdb"))[0]  # Take only the first model. Trivial for X-ray, only 1st for NMR
            io = PDBIO()
            io.set_structure(structure)
            io.save(os.path.join(os.path.join(path, label, label + ".pdb")), DomSelect())


            ## 4. Select ligands and save them separately
            pdb_structure = parser.get_structure("st", os.path.join(path, label, pdb + ".pdb"))[0]  # Take only the first model. Trivial for X-ray, only 1st for NMR
            ligands = [i for i in pdb_structure.get_residues() if i.get_resname() == intlig]  # Repeated entities of the same ligand (e.g 'YDJ repeated twice')

            class LigSelect(Select):
                def accept_residue(self, residue):
                    if residue == ligand:
                        return 1
                    else:
                        return 0

            def SASA(prot, lig):
    
                # compute ligand SASA
                lig_h = Chem.AddHs(lig, addCoords=True)

                # Get Van der Waals radii (angstrom)
                ptable = Chem.GetPeriodicTable()
                radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in lig_h.GetAtoms()]

                # Compute solvent accessible surface area
                lig_sasa = rdFreeSASA.CalcSASA(lig_h, radii)
                
                # Join protein & ligand
                comp = Chem.CombineMols(prot, lig)
                comp_h = Chem.AddHs(comp, addCoords=True)

                # Get Van der Waals radii (angstrom)
                ptable = Chem.GetPeriodicTable()
                radii = [ptable.GetRvdw(atom.GetAtomicNum()) for atom in comp_h.GetAtoms()]

                # Compute solvent accessible surface area
                comp_sasa = rdFreeSASA.CalcSASA(comp_h, radii)
                
                comp_lig = Chem.GetMolFrags(comp_h, asMols=True,  sanitizeFrags=True)
                comp_lig = [i for i in comp_lig if MolWt(i) == MolWt(lig_h)]
                # Avoid including smaller ligands than expected (e.g. 7PE in 5EHY)
                at = [i.GetAtomicNum() for i in lig.GetAtoms()]
                if len(comp_lig) != 1 or Counter(at)[6] < 5:
                    return np.nan, np.nan
                else:
                    comp_lig = comp_lig[0]

                lig_sasa_free = 0
                for a in lig_h.GetAtoms():
                    lig_sasa_free += float(a.GetProp("SASA"))

                lig_sasa_bound = 0
                for a in comp_lig.GetAtoms():
                    lig_sasa_bound += float(a.GetProp("SASA"))

                return round(lig_sasa_free, 3), round(lig_sasa_bound, 3)



            for count, ligand in enumerate(ligands):

                try:

                    # Save ligand
                    io = PDBIO()
                    io.set_structure(pdb_structure)
                    io.save(os.path.join(os.path.join(path, label, intlig + "_" + str(count) + ".pdb")), LigSelect())
                    #command = 'obabel ' + os.path.join(path, label, intlig + "_" + str(count) + ".pdb") + ' -O ' + os.path.join(path, label, intlig + "_" + str(count) + ".sdf")
                    #os.system(command)
                    #os.remove(os.path.join(path, label, intlig + "_" + str(count) + ".pdb"))

                    # Save centroid
                    res = [i for i in pdb_structure.get_residues() if i == ligand][0]
                    ligatoms = [at.coord for at in res.get_atoms()]
                    x = np.mean(np.array(ligatoms)[:,0])
                    y = np.mean(np.array(ligatoms)[:,1])
                    z = np.mean(np.array(ligatoms)[:,2])
                    center = np.array([x, y, z], dtype=np.float32)

                    x, y, z = str(round(x, 3)), str(round(y, 3)), str(round(z, 3))
                    ctr = " "*(8-len(x)) + x + " "*(8-len(y)) + y + " "*(8-len(z)) + z
                    text = """HEADER\nHETATM    1   C  CTR A   1    """ + ctr + """  1.00  1.00           C\nEND"""

                    with open(os.path.join(path, label, intlig + "_CTR_" + str(count) + ".pdb"), "w") as f:
                        f.write(text)


                    # Get SASA
                    prot = Chem.MolFromPDBFile(os.path.join(path, label, label + ".pdb"))
                    lig = Chem.MolFromPDBFile(os.path.join(path, label, intlig + "_" + str(count) + ".pdb"))

                    lig_sasa_free, lig_sasa_bound = SASA(prot, lig)
                    acc = round(lig_sasa_bound/lig_sasa_free, 3)

                    outfile.write("\t".join([domain, pdb, chain, intlig, str(count), str(lig_sasa_free), str(lig_sasa_bound), str(acc)]) + "\n")

                except:

                    outfile.write("\t".join([domain, pdb, chain, intlig, str(count), str('failed'), str('failed'), str('failed')]) + "\n")


        os.chdir(path)   
        tar = tarfile.open(label + ".tar.gz", "w:gz")
        tar.add(label)
        shutil.rmtree(label)
        tar.close()