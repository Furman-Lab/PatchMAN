## PatchMAN protocol for blind peptide-protein docking

<img align="left" width="500" height="144" src="https://raw.githubusercontent.com/Alisa-Kh/PatchMAN/master/img/PatchMAN_small.PNG">

<br /><br /><br /><br /><br />

### Description

PatchMAN (Patch-Motif AligNments) maps the receptor surface for local structural motif matches in structures of protein monomers and interfaces, to extract complementary fragments and derive templates for peptide-protein docking.

The protocol consists of 4 consecutive steps: (1) Definition of surface patches on the receptor; (2) Identification of structural motif matches in protein structures, and an interacting fragment that can be used as template for the bound peptide; (3) Generation of the peptide-protein complex template structure, by superimposing the interacting peptide back onto the receptor, and (4) Replacing side chains according to the peptide sequence (threading), refinement and scoring of the model.

---

### Installation

#### Downloading software and data
1. Obtain [Rosetta](https://els2.comotion.uw.edu/product/rosetta) and [PyRosetta](https://els2.comotion.uw.edu/product/pyrosetta) licenses.
2. Register for [MASTER v1.6](https://grigoryanlab.org/index.php?sec=get&soft=MASTER) and obtain download URL.
3. Create a .env file to contain login information: ```cp sample.env .env```
4. Edit `.env` file with the Rosetta and PyRosetta usernames and passwords and MASTER's URL.
5. Run `bash download_data_and_software.sh` that will take care of downloading all the required files and extracting the database for MASTER search into the `databases/masterDB` and `databases/master_clean` directory. Downloading Rosetta, PyRosetta and the MASTER database can take significant amount of time, depending on your network. The script downloads the versions of the softwares that are used in the paper.

:exclamation: Note: Running MPI in Singularity containers require that the version of the hose and container MPI match. The script automatically detects the version of host OpenMPI (required to speed up FlexPepDock runs) and downloads it. If the container is not built on the host computer that will run it, the variable OMPI_VERSION might need to be manually modified.


#### Option 1: Installation with Singularity containers (recommended)

Three singularity definition files are provided for compiling MASTER, Rosetta and PyRosetta together with all the required python scripts. 
How you need to build singularity images might be system dependent, as it requires sudo. For example, you might need to run virtual machines or other similar systems. The images can be built with:
```
sudo singularity build rosetta.sif rosetta.def # compiling Rosetta can take significant amount of time
sudo singularity build python.sif python.def
sudo singularity build master.sif master.def # this also takes care of patching PatchMAN
```
#### Option 2: Installation without containers

- To install python packages and [PyRosetta](https://www.pyrosetta.org/downloads), create a virtual environment and activate it. One example:
    ```
    virtualenv patchman --python=python3
    source patchman/bin/activate  # or activate.csh, based on your shell type
    pip3 install -r requirements.txt # install required packages
    python3 --version # get your Python3 version
    ```
    Download PyRosetta wheel that matches your python3 version and your OS with the link, username and password sent [after registration](https://els2.comotion.uw.edu/product/pyrosetta).
    ```
    pip3 install <downloaded pyrosetta wheel> 
    ```
    To get support on PyRosetta and its installation, visit: https://www.pyrosetta.org/downloads
- Install Rosetta from the downloaded rosetta.tar.gz with mpi support. You can change the number of used cores (-j argument) according to your system
    ```
    cd containers/
    mkdir rosetta
    tar -xzf rosetta.tar.gz -C rosetta --strip-components=1
    cd rosetta/main/source 
    python scons.py -j 4 extras=mpi,serialization mode=release bin/FlexPepDocking.mpiserialization.linuxgccrelease bin/cluster.mpiserialization.linuxgccrelease
    ```
    For further support on Rosetta installation, please refer to the documentation.
- Set up MASTER
  - The downloaded source code needs a slight modification for running PatchMAN . This can be done in 2 different ways:
    -  Programatically (using the patch file in the `bin/` directory):
        ```
        patch -l master/src/Match.cpp bin/master.patch
        ```
    -  Manually: Go to line 107 in the `Match.cpp` file and add the following code (before `return os;`):

        ```
        double *T=((Match*)(&m))->getTranslation();
        double **R=((Match*)(&m))->getRotation();

        os << " T: " << T[0]    << " " << T[1]    << " " << T[2]    << " ";
        os << " U: " << R[0][0] << " " << R[0][1] << " " << R[0][2] << " "
                     << R[1][0] << " " << R[1][1] << " " << R[1][2] << " "
                     << R[2][0] << " " << R[2][1] << " " << R[2][2] << " ===" ;
        ```
  - Follow the instruction in the INSTALL file to compile MASTER.
---

### Quick start

PatchMAN can be run with:

`bash PatchMAN_protocol.sh <arguments> RECEPTOR PEPTIDE`

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
*bioRxiv 2021.09.02.458699; doi:https://doi.org/10.1101/2021.09.02.458699*  

Please also cite the following papers:

**Rapid Search for Tertiary Fragments Reveals Protein Sequence-Structure Relationships**  
Zhou J., Grigoryan G.  
*Protein Science, 24(4): 508-524, 2015.*  

**Sub-angstrom modeling of complexes between flexible peptides and globular proteins**  
Raveh B, London N, Schueler-Furman O. (2010).  
*Proteins 78:2029-40.*  

**PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta**  
Chaudhury S, Lyskov S, Gray JJ.  
*Bioinformatics. 2010;26(5):689-691. doi:10.1093/bioinformatics/btq007*  
