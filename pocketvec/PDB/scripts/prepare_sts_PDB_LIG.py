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



def get_xml(pdb_id, path_to_xml='.', force_rerun=False):

    """
    Inspired in ssbio, def download_sifts_xml

    """
    
    filename = '{}.xml.gz'.format(pdb_id.lower())
    url = os.path.join('http://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/', filename)
    outfile = os.path.join(path_to_xml, pdb_id.lower() + '.xml.gz')
    out_xml = os.path.join(path_to_xml, pdb_id.lower() + '_' + str(task_id) + '.xml')
    
    # Download xml.gz
    if force_rerun is True or pdb_id.lower() + '.xml.gz' not in set(os.listdir(path_to_xml)):
        urllib.request.urlretrieve(url, outfile)
        
    # Uncompress
    with gzip.open(outfile, 'rb') as f:
        with open(out_xml, "wb") as out:
            file_content = f.read()
            out.write(file_content)

    time.sleep(1)
            
    return out_xml



def map_pdbres_to_uniprot(pdb_res, chain_id, sifts_file):
    
    """
    Given a pdb residue, a chain id and a sifts_file, map it with the corresponding
    uniprot residue   
    
    """

    parser = etree.XMLParser(ns_clean=True)
    tree = etree.parse(sifts_file, parser)
    root = tree.getroot()

    my_pdb_resnum = None

    # TODO: "Engineered_Mutation is also a possible annotation, need to figure out what to do with that
    my_pdb_annotation = False

    # Find the right chain (entities in the xml doc)
    ent = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}entity'


    for chain in root.findall(ent):
        if True: #chain.attrib['entityId'] == chain_id: # this is not right (e.g. 1a0h, chain E in pdb is chain D in uniprot! check xml)
            # Find the "crossRefDb" tag that has the attributes dbSource="PDB" and  dbResNum="pdb_res"
            # Then match it to the crossRefDb dbResNum that has the attribute dbSource="PDBresnum"
            # Check if uniprot + resnum even exists in the sifts file (it won't if the pdb doesn't contain the residue)
            ures = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="PDB"][@dbResNum="%s"][@dbChainId="%s"]' % (pdb_res, chain_id)
            my_pdb_residue = chain.findall(ures)
            
            if len(my_pdb_residue) == 1:
                # Get crossRefDb dbSource="PDB"
                parent = my_pdb_residue[0].getparent()
                pres = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}crossRefDb[@dbSource="UniProt"]'
                my_pdb_residue = parent.findall(pres)
                if len(my_pdb_residue) == 0 or my_pdb_residue[0].attrib['dbResNum'] == 'null':
                    pass
                else:
                    my_pdb_resnum = int(my_pdb_residue[0].attrib['dbResNum'])


                    # Get <residueDetail dbSource="PDBe" property="Annotation">
                    # Will be Not_Observed if it is not seen in the PDB
                    anno = './/{http://www.ebi.ac.uk/pdbe/docs/sifts/eFamily.xsd}residueDetail[@dbSource="PDBe"][@property="Annotation"]'
                    my_pdb_annotation = parent.findall(anno)
                    if len(my_pdb_annotation) == 1:
                        my_pdb_annotation = my_pdb_annotation[0].text
                        if my_pdb_annotation == 'Not_Observed':
                            my_pdb_annotation = False
                    else:
                        my_pdb_annotation = True

                    return pdb_res, my_pdb_resnum, my_pdb_annotation
            else:
                pass

    return pdb_res, my_pdb_resnum, my_pdb_annotation





with open(path_to_summary, "w") as outfile:


    for info in infos:


        try:


            pdb = info[0]
            int_lig = info[1]
            lig = info[2]
            het = info[3]
            chain_id = info[4]
            uniprot = info[5]

            path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/structures", pdb[1:3], pdb + "_" + chain_id + "_" + int_lig)
            if os.path.exists(path) is False: os.makedirs(path)
            os.chdir(os.path.join(path))

            ### 1. Create directory & download pdb structure
            with gzip.open('/aloy/home/acomajuncosa/programs/localpdb/mirror/pdb/' + pdb[1:3].lower() + '/pdb' + pdb.lower() + '.ent.gz', 'rb') as f_in:
                with open(os.path.join(path, pdb + ".pdb"), 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                    
            ### Check ligands to remove: INT_LIG, LIG but not HET (nonstd residues)
            ligands_to_remove = ['W']
            ligands_to_remove.append("H_" + " "*(3-len(int_lig)) + int_lig)
            ligands_to_remove.append("H_" + int_lig + " "*(3-len(int_lig)))

            if type(lig) == list:
                for l in lig:
                    ligands_to_remove.append("H_" + " "*(3-len(l)) + l)
                    ligands_to_remove.append("H_" + l + " "*(3-len(l)))
                ligands_to_remove = set(ligands_to_remove)



            ### 2. Take only the first model
            parser = PDBParser()
            path_in = os.path.join(path, pdb + ".pdb")
            path_out = path
            structure = parser.get_structure("st", path_in) # Take only the first model. Trivial for X-ray, only 1st for NMR

            if len(structure) > 1:
                structure = structure[0]

            io = PDBIO()
            io.set_structure(structure)
            io.save(os.path.join(path_out, pdb + "_model.pdb"))


            ### 3. Remove water molecules and ligands, but not nonstd residues, and select the studied chain
            parser = PDBParser()
            path_in = os.path.join(path, pdb + "_model.pdb")
            path_out = path
            structure = parser.get_structure("st", path_in)

            class remove_ligs(Select):
                def accept_residue(self, residue):
                    if residue.get_id()[0] not in ligands_to_remove and residue.get_parent().id == chain_id:
                        return 1
                    else:
                        return 0

            io = PDBIO()
            io.set_structure(structure)
            io.save(os.path.join(path, pdb + "_" + chain_id + "_st.pdb"), remove_ligs())


            ### 4. Get the number of interesting ligands in the chain. If none, do not continue
            interesting_ligands = [i for i in structure.get_residues() if i.get_resname() == int_lig and i.get_parent().id == chain_id]

            if len(interesting_ligands) > 0:


                ### 5. Remove hydrogens
                command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + os.path.join(path_out , pdb + "_" + chain_id + '_st.pdb') + ' -o ' + os.path.join(path_out , pdb + "_" + chain_id + '_hydrogens.pdb') + ' --force_save --non_interactive rem_hydrogen --remove Yes'
                o = os.popen(command).read()
                sys.stderr.write(o + "\n\n")
                sys.stderr.flush()

                ### 6. Select occupancies
                command = 'python /aloy/home/acomajuncosa/programs/structureChecking/bin/check_structure -i ' + os.path.join(path_out , pdb + "_" + chain_id + '_hydrogens.pdb') + ' -o ' + os.path.join(path_out , pdb + "_" + chain_id + '_altloc.pdb') + ' --force_save --non_interactive altloc --select occupancy'
                o = os.popen(command).read()
                sys.stderr.write(o + "\n\n")
                sys.stderr.flush()

                shutil.copyfile(os.path.join(path_out , pdb + "_" + chain_id + '_altloc.pdb'), os.path.join(path_out , pdb + "_" + chain_id + '_' + int_lig + '.pdb'))


                try:

                    ### 7. Check coverage

                    ### 7.1. Get xml
                    path_to_xml = "/aloy/home/acomajuncosa/PocketVec_v2/HT/data/PDB/coverage/xml"
                    path_to_sifts = get_xml(pdb, path_to_xml, force_rerun=False)

                    ### 7.2 Get uniprot sequence
                    baseUrl = "http://www.uniprot.org/uniprot/"
                    currentUrl = baseUrl + uniprot + ".fasta"
                    response = r.post(currentUrl)
                    cData = ''.join(response.text)
                    Seq = StringIO(cData)
                    pSeq = list(SeqIO.parse(Seq,'fasta'))
                    prot_sequence = str(pSeq[0].seq)


                    ### 7.3 Get coverage residue by residue


                    parser = PDBParser()
                    structure = parser.get_structure("st", os.path.join(path_out , pdb + '_' + chain_id + "_" + int_lig + '.pdb'))
                    cov = {}
                    for i in structure.get_residues():
                        numb = str(i.get_id()[1])
                        a, b, c = map_pdbres_to_uniprot(numb, chain_id, path_to_sifts)
                        cov[a] = [b, c]

                    os.remove(path_to_sifts)

                    numb_covered = Counter([cov[i][1] for i in cov])[True]
                    total_seq = len(prot_sequence)
                    coverage = round(numb_covered/total_seq, 3)

                except:

                    coverage = -1  ## Coverage calculation failed

                if coverage >= 0.8:

                    ### 8. Save ligand independently and calculate SASA (free/bound)

                    parser = PDBParser()
                    structure = parser.get_structure("st", os.path.join(path, pdb + "_model.pdb"))

                    for c, interesting_ligand in enumerate(interesting_ligands):
                        
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
                        
                        # Get centroid
                        get_centroid(os.path.join(path, int_lig + "_" + str(c) + ".pdb"), os.path.join(path, int_lig + "_" + str(c) + "_centroid.pdb"))
                        command = 'obabel ' + os.path.join(path, int_lig + "_" + str(c) + "_centroid.pdb") + " -O " + os.path.join(path, int_lig + "_" + str(c) + "_centroid.sd")
                        os.system(command)

                        try:
                        
                            # Get SASA
                            prot = Chem.MolFromPDBFile(os.path.join(path, pdb + "_" + chain_id + "_" + int_lig + ".pdb"))
                            lig = Chem.MolFromPDBFile(os.path.join(path, int_lig + "_" + str(c) + ".pdb"))
                            if lig.GetNumHeavyAtoms() == Chem.MolFromInchi(pdbcode_to_inchi[int_lig]).GetNumHeavyAtoms():
                                lig_sasa_free, lig_sasa_bound = SASA(prot, lig)
                            else:
                                lig_sasa_free, lig_sasa_bound = np.nan, np.nan
                            acc = round(lig_sasa_bound / lig_sasa_free, 3)

                            outfile.write("\t".join([pdb, chain_id, uniprot, int_lig, str(coverage), str(c), str(acc), 'all good']) + "\n")  # Everything fine (caution, acc can be nan here)

                        except:

                            outfile.write("\t".join([pdb, chain_id, uniprot, int_lig, str(coverage), str(c), "error", 'SASA calculation failed']) + "\n")  # SASA calculation failed

                else:

                    outfile.write("\t".join([pdb, chain_id, uniprot, int_lig, str(coverage), "--", "--", 'coverage < 0.8']) + "\n")  # Coverage < 0.8


            else:

                outfile.write("\t".join([pdb, chain_id, uniprot, int_lig, "--", "--", "--", 'ligand not in chain']) + "\n")  # Ligand not in chain


        except:

            outfile.write("\t".join([pdb, chain_id, uniprot, int_lig, "--", "--", "--", "first steps error"]) + "\n")  # Failed in the starting steps



        path = os.path.join("/aloy/home/acomajuncosa/MurD/GitHub/structures", pdb[1:3])
        os.chdir(path)   
        tar = tarfile.open(pdb + "_" + chain_id + "_" + int_lig + ".tar.gz", "w:gz")
        tar.add(pdb + "_" + chain_id + "_" + int_lig)
        shutil.rmtree(pdb + "_" + chain_id + "_" + int_lig)
        tar.close()
