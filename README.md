## PatchMAN protocol for blind peptide-protein docking

<img align="left" width="500" height="144" src="https://raw.githubusercontent.com/Alisa-Kh/PatchMAN/master/img/PatchMAN_small.PNG"> 

<br /><br /><br /><br /><br />

### Description

PatchMAN (Patch-Motif AligNments) maps the receptor surface for local structural motif matches in structures of protein monomers and interfaces, to extract complementary fragments and derive templates for peptide-protein docking.

The protocol consists of 4 consecutive steps: (1) Definition of surface patches on the receptor; (2) Identification of structural motif matches in protein structures, and an interacting fragment that can be used as template for the bound peptide; (3) Generation of the peptide-protein complex template structure, by superimposing the interacting peptide back onto the receptor, and (4) Replacing side chains according to the peptide sequence (threading), refinement and scoring of the model.

An example protocol, set up using Singularity containers with Slurm workload manager can be found in PatchMAN_protocol.sh. 

All the softwares and data can be downloaded with `download_data_and_software.sh`. The usernames and passwords needed to download Rosetta and PyRrosetta are read from a .env file that needs to be generated based on `sample.env`, after obtaining licenses for both of them. The version of MPI on the host system is extracted automatically

---  

### Singularity installation  

Master search, extraction of peptides based on the structural motif search results and Rosetta FlexPepDock runs are best to be run distributed, with MPI. The protocol script is setup to use OpenMPI with Slurm,

The provided Singularity files need to be build with the following commands, after downloading the needed files with `download_data_and_software.sh`:
```
singularity build rosetta.sif rosetta.def
singularity build python.sif python.def
singularity build master.sif master.def
```

After builds are complete (Rosetta and PyRosetta can take some time) running PatchMAN with 
`bash PatchMAN_protocol.sh <arguments> RECEPTOR PEPTIDE`
where RECEPTOR is a PDB file and PEPTIDE is the peptide sequence to be docked. The peptide can contain post-translational modifications, denoted by Rosetta standards, e.g. [SER:phosphorylated]. The available PTMs can be listed with `singularity run rosetta.sif ls /rosetta/rosetta_src_2019.14.60699_bundle/main/database/chemical/residue_type_sets/fa_standard/patches`

---  

### Command line installation

#### Software prequisites

To run PatchMAN the following prequisites should be downloaded and installed:

1. [Python(3.5)](https://www.python.org/downloads/source/)
2. [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download)
3. [MASTER v1.6](https://grigoryanlab.org/master/)
4. [Rosetta](https://www.rosettacommons.org/software/license-and-download)


#### Installation

- Download and install [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download) and [Rosetta](https://www.rosettacommons.org/software/license-and-download)
- Set up MASTER
  - Download the source code of [MASTER v1.6](https://grigoryanlab.org/index.php?sec=get&soft=MASTER)
  - In the Match.cpp file:
Go to line 107 and add the following code (before ```return os;```):
```
    double *T=((Match*)(&m))->getTranslation();
    double **R=((Match*)(&m))->getRotation();

    os << " T: " << T[0]    << " " << T[1]    << " " << T[2]    << " ";
    os << " U: " << R[0][0] << " " << R[0][1] << " " << R[0][2] << " "
                 << R[1][0] << " " << R[1][1] << " " << R[1][2] << " "
                 << R[2][0] << " " << R[2][1] << " " << R[2][2] << " ===" ;
```
  - Follow the instruction in the INSTALL file to compile MASTER 
- Download [MASTER database](https://grigoryanlab.org/master/#database) for template search

#### Running PatchMAN

A PatchMAN demo run can be found in the *example_run* folder
