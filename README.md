## PatchMAN protocol for blind peptide-protein docking

### Description

PatchMAN (Patch-Motif AligNments), a novel approach that maps the receptor surface for local structural motif matches in structures of protein monomers and interfaces, to extract complementary fragments and derive templates for peptide-protein docking.

The protocol consists of 4 consecutive steps: (1) Definition of surface patches on the receptor; (2) Identification of structural motif matches in protein structures, and an interacting fragment that can be used as template for the bound peptide; (3) Generation of the peptide-protein complex template structure, by superimposing the interacting peptide back onto the receptor, and (4) Replacing side chains according to the peptide sequence (threading), refinement and scoring of the model.


### Software prequisites

To run PatchMAN the following prequisites should be downloaded and installed:

1. [Python(3.5)](https://www.python.org/downloads/source/)
2. [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download)
3. [MASTER](https://grigoryanlab.org/master/)
4. [Rosetta](https://www.rosettacommons.org/software/license-and-download)

### Running PatchMAN

PatchMAN demo run can be found in the *example_run* folder
