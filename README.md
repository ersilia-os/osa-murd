# New ligands for the MurD Open Source Antibiotics challenge
Open Source Antibiotics MurD ligands based on generative models.

This repository contains scripts and data used to produce MurD ligand candidates for Open Source Antibiotics (OSA), as described in issue [#69](https://github.com/opensourceantibiotics/murligase/issues/69). Most information realted to MurD can be found in the corresponding [murligase](https://github.com/opensourceantibiotics/murligase) OSA repository.

## The challenge
The target is the bacterial enzyme [MurD](https://github.com/opensourceantibiotics/murligase/wiki). We have been asked to take existing crystal structures of MurD and generate novel small molecules tha tbind to an allosteric pocket.

## Available data

### Open Source Antibiotics
* Structures of the four fragments (349, 373, 374, 378) bound to the **allosteric** MurD site of interest: [Fragalysis](https://fragalysis.diamond.ac.uk/viewer/react/preview/target/MURD), [Wiki](https://github.com/opensourceantibiotics/murligase/wiki/Initial-MurD-Hits), [PDB files](https://github.com/opensourceantibiotics/murligase/tree/master/docs/pdbs_forNGL/MurD). Structures correspond to *Streptococcus agalactiae* MurD.
* An *Escherichia coli* MurD structure (PDB [3UAG](https://www.rcsb.org/structure/3UAG)) without fragments.
* The *Spreptococcus agalactiae* MurD structure (PDB [3LK7](https://www.rcsb.org/structure/3lk7)).
* Experimental [inhibition measurements](https://github.com/opensourceantibiotics/murligase/wiki/MurD-Round-1) of some fragments.
* [Atomwise predictions](https://github.com/opensourceantibiotics/murligase/tree/master/Atomwise) of binding to the 3UAG structure (see OSA issue [#33](https://github.com/opensourceantibiotics/murligase/issues/33)).
* Crystal structures of some of Atomwise predictions, but [bound to a related protein, MurE](https://github.com/opensourceantibiotics/murligase/wiki/XChem-EcMurE-Atomwise-library-1).

### From other sources

## Our approach
