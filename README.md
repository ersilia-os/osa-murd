# New ligands for the MurD Open Source Antibiotics challenge
Open Source Antibiotics MurD ligands based on generative models.

This repository contains scripts and data used to produce MurD ligand candidates for Open Source Antibiotics (OSA), as described in issue [#69](https://github.com/opensourceantibiotics/murligase/issues/69). Most information realted to MurD can be found in the corresponding [murligase](https://github.com/opensourceantibiotics/murligase) OSA repository.

## The challenge
The target is the bacterial enzyme [MurD](https://github.com/opensourceantibiotics/murligase/wiki). We have been asked to take existing crystal structures of MurD and generate novel small molecules that bind to an allosteric pocket.

## Available data

### Open Source Antibiotics
* Structures of the four fragments (349, 373, 374, 378) bound to the **allosteric** MurD site of interest: [Fragalysis](https://fragalysis.diamond.ac.uk/viewer/react/preview/target/MURD), [Wiki](https://github.com/opensourceantibiotics/murligase/wiki/Initial-MurD-Hits), [PDB files](https://github.com/opensourceantibiotics/murligase/tree/master/docs/pdbs_forNGL/MurD). Structures correspond to *Streptococcus agalactiae* MurD.
* An *Escherichia coli* MurD structure (PDB [3UAG](https://www.rcsb.org/structure/3UAG)) without fragments.
* The *Spreptococcus agalactiae* MurD structure (PDB [3LK7](https://www.rcsb.org/structure/3lk7)).
* Experimental [inhibition measurements](https://github.com/opensourceantibiotics/murligase/wiki/MurD-Round-1) of some fragments.
* [Atomwise predictions](https://github.com/opensourceantibiotics/murligase/tree/master/Atomwise) of binding to the 3UAG structure (see OSA issue [#33](https://github.com/opensourceantibiotics/murligase/issues/33)).
* Crystal structures of some of Atomwise predictions, but [bound to a related protein, MurE](https://github.com/opensourceantibiotics/murligase/wiki/XChem-EcMurE-Atomwise-library-1).

## Results

Our suggested molecules can be found in the Open Source Antibiotics GitHub page, issue [#76](https://github.com/opensourceantibiotics/murligase/issues/79). Additionally, we have generated a [Molecular Cloud](https://github.com/whitead/molcloud) made up by our candidate molecules:

![molecular cloud](https://user-images.githubusercontent.com/80755454/177553868-bb44f903-73c1-42b8-a6c1-81b887c3ff93.jpeg)


## About us

This work is the result of a collaboration between the [Structural Bioinformatics and Network Biology Laboratory](https://sbnb.irbbarcelona.org), at IRB Barcelona and the [Ersilia Open Source Initiative](https://ersilia.io). 
