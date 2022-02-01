## PatchMAN protocol for blind peptide-protein docking

<img align="left" width="500" height="144" src="https://raw.githubusercontent.com/Alisa-Kh/PatchMAN/master/img/PatchMAN_small.PNG"> 

<br /><br /><br /><br /><br />

### Description

PatchMAN (Patch-Motif AligNments) maps the receptor surface for local structural motif matches in structures of protein monomers and interfaces, to extract complementary fragments and derive templates for peptide-protein docking.

The protocol consists of 4 consecutive steps: (1) Definition of surface patches on the receptor; (2) Identification of structural motif matches in protein structures, and an interacting fragment that can be used as template for the bound peptide; (3) Generation of the peptide-protein complex template structure, by superimposing the interacting peptide back onto the receptor, and (4) Replacing side chains according to the peptide sequence (threading), refinement and scoring of the model.

---  

### Installation  

#### Downloading software and data
1. Obtain [Rosetta](https://www.rosettacommons.org/software/license-and-download) and [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download) licenses.
2. Create a .env file to contain login information: ```cp sample.env .env```
3. Edit `.env` file with the Rosetta and PyRosetta usernames and passwords.
4. Run `bash download_data_and_software.sh` that will take care of downloading all the required files and extracting the database for MASTER search into the `masterDB` directory. Downloading Rosetta, PyRosetta and the MASTER database can take significant amount of time, depending on your network. The script downloads the versions of the softwares that are used in the paper.

:exclamation: Note: Running MPI in Singularity containers require that the version of the hose and container MPI match. The script automatically detects the version of host OpenMPI (required to speed up FlexPepDock runs) and downloads it. If the container is not built on the host computer that will run it, the variable OMPI_VERSION might need to be manually modified.


#### Installation with Singularity containers

Three singularity definition files are provided for compiling MASTER, Rosetta and PyRosetta together with all the required python scripts. The images can be built with

```
sudo singularity build rosetta.sif rosetta.def # compiling Rosetta can take significant amount of time
sudo singularity build python.sif python.def
sudo singularity build master.sif master.def # this also takes care of patching PatchMAN
```

#### Installation without containers

- The required python packages can be installed with `pip install -r requirements.txt`
- Install [PyRosetta](https://www.pyrosetta.org/downloads/legacy-pyrosetta3-download) and [Rosetta](https://www.rosettacommons.org/software/license-and-download) from the downloaded rosetta.tar.gz and wheel files.
- Set up MASTER
  - Download the source code of [MASTER v1.6](https://grigoryanlab.org/index.php?sec=get&soft=MASTER). The source code needs a slight modification for running PatchMAN. This can be done in 2 different ways:
    -  Programatically: use the patch file in the `bin/` directory and issue:
    
        ```patch -l path_to_master/src/Match.cpp bin/master.patch```
    
    -  Manually: Go to line 107 in the `Match.cpp` file and add the following code (before ```return os;```):
    
        ```
        double *T=((Match*)(&m))->getTranslation();
        double **R=((Match*)(&m))->getRotation();

        os << " T: " << T[0]    << " " << T[1]    << " " << T[2]    << " ";
        os << " U: " << R[0][0] << " " << R[0][1] << " " << R[0][2] << " "
                     << R[1][0] << " " << R[1][1] << " " << R[1][2] << " "
                     << R[2][0] << " " << R[2][1] << " " << R[2][2] << " ===" ;
        ```
    
  - Follow the instruction in the INSTALL file to compile MASTER 

---

### Quick start

PatchMAN can be run with: 

```bash PatchMAN_protocol.sh <arguments> RECEPTOR PEPTIDE```  

where RECEPTOR is a PDB file and PEPTIDE is the peptide sequence to be docked.   
The peptide can contain post-translational modifications, denoted by Rosetta standards, e.g. `[SER:phosphorylated]`. The available PTMs can be listed with   

```singularity run rosetta.sif ls /rosetta/main/database/chemical/residue_type_sets/fa_standard/patches```  

:exclamation: Note that the protocol script is set up to use Singularity containers. If you compiled Rosetta, PyRosetta or MASTER without containers, you will need to edit the `$ROSETTA`, `$PYTHON` and `$MASTER` environmental variables accordingly.

:exclamation: Note that the protocol script is set up to use Slurm job scheduler. Using an other type of scheduler needs editing of the `PatchMAN_protocol.sh` file and `.sh` files in the `bin/` directory. Unfortunately, we cannot help with that.


#### Test run
A test run of PatchMAN can be performed on the 1ssh.pdb in the `test/` directory. Turning off receptor backbone minimization for testing purposes decreases runtime: 

```
cd test/
../PatchMAN_protocol.sh -m false 1ssh.pdb EGPPPAMPARPT
```

### Running parameters
```
-m minimize receptor backbone (true or false; default: true)
-t number of refinement runs for FlexPepDock (default: 1)
-s masked residues as a PDB structure (default: None) # Note: this feature is experimental

-w working directory (default: current directory)
-o output directory for the ligands (default: working directory)
-g log file (Default is stdout)
-e error log file (Default is stderr)
-n job name
-v verbose (default: false)
```

---
### Citing PatchMAN

**PatchMAN docking: Modeling peptide-protein interactions in the context of the receptor surface**  
Alisa Khramushin, Tomer Tsaban, Julia Varga, Orly Avraham, Ora Schueler-Furman  
*bioRxiv 2021.09.02.458699; doi: https://doi.org/10.1101/2021.09.02.458699*  
