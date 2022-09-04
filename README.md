## PatchMAN protocol for blind peptide-protein docking

<img align="left" width="500" height="144" src="https://raw.githubusercontent.com/Alisa-Kh/PatchMAN/master/img/PatchMAN_small.PNG"> 

<br /><br /><br /><br /><br />


### Description

PatchMAN (Patch-Motif AligNments) maps the receptor surface for local structural motif matches in structures of protein monomers and interfaces, to extract complementary fragments and derive templates for peptide-protein docking.

The protocol consists of 4 consecutive steps: (1) Definition of surface patches on the receptor; (2) Identification of structural motif matches in protein structures, and an interacting fragment that can be used as template for the bound peptide; (3) Generation of the peptide-protein complex template structure, by superimposing the interacting peptide back onto the receptor, and (4) Replacing side chains according to the peptide sequence (threading), refinement and scoring of the model.


### Software prequisites

To run PatchMAN the following prequisites should be downloaded and installed:

1. [Python(3.5)](https://www.python.org/downloads/source/)
2. [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download)
3. [MASTER v1.6](https://grigoryanlab.org/master/)
4. [Rosetta](https://www.rosettacommons.org/software/license-and-download)

### Installation

- Download and install [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download) and [Rosetta](https://www.rosettacommons.org/software/license-and-download)
- Set up MASTER
  - Download the source code of [MASTER v1.6](https://grigoryanlab.org/index.php?sec=get&soft=MASTER)
  - In the Match.cpp file:
Go to line 107 and add the following code (before ```return os;```):

    ```double *T=((Match*)(&m))->getTranslation();
    double **R=((Match*)(&m))->getRotation();

    os << " T: " << T[0]    << " " << T[1]    << " " << T[2]    << " ";
    os << " U: " << R[0][0] << " " << R[0][1] << " " << R[0][2] << " "
                 << R[1][0] << " " << R[1][1] << " " << R[1][2] << " "
                 << R[2][0] << " " << R[2][1] << " " << R[2][2] << " ===" ;```
  - Follow the instruction in the INSTALL file to compile MASTER 
- Download [MASTER database](https://grigoryanlab.org/master/#database) for template search

### Running PatchMAN

A PatchMAN demo run can be found in the *example_run* folder

### Citing PatchMAN

PatchMAN docking: Modeling peptide-protein interactions in the context of the receptor surface
Alisa Khramushin, Tomer Tsaban, Julia Varga, Orly Avraham, Ora Schueler-Furman
bioRxiv 2021.09.02.458699; doi: https://doi.org/10.1101/2021.09.02.458699

