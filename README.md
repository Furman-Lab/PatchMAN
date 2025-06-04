## PatchMAN protocol for blind peptide-protein docking

<img align="left" width="500" height="144" src="https://raw.githubusercontent.com/Alisa-Kh/PatchMAN/master/img/PatchMAN_small.PNG">

<br /><br /><br /><br /><br />

### Description

PatchMAN (Patch-Motif AligNments) maps the receptor surface for local structural motif matches in structures of protein monomers and interfaces, to extract complementary fragments and derive templates for peptide-protein docking.

The protocol consists of 4 consecutive steps: (1) Definition of surface patches on the receptor; (2) Identification of structural motif matches in protein structures, and an interacting fragment that can be used as template for the bound peptide; (3) Generation of the peptide-protein complex template structure, by superimposing the interacting peptide back onto the receptor, and (4) Replacing side chains according to the peptide sequence (threading), refinement and scoring of the model.

---

### Installation

#### Downloading software and data
1. Obtain the files in this repo: 
   `git clone --branch python https://github.cs.huji.ac.il/fora-lab/PatchMAN-forServer`  
2. Register for [MASTER v1.6](https://grigoryanlab.org/index.php?sec=download&soft=MASTER), go to the download page and copy the link of the `Source code` of MASTER v1.6 and add its URL to create_env.sh (to define $MASTER_URL).
3. Run `bash download_data_and_software.sh`, that will take care of downloading all the required files and extracting the databases to the `databases/` directory. This usually takes 5-10mins.

#### Installation with conda/mamba environment (recommended)
To create a new environment, just simply run the following command:
```
bash create_env.sh
```
:exclamation: the script uses micromamba for Linux. If you use another OS, change the `MICROMAMBA_URL` accordingly and might need to modify a few commands also. See:
- MacOS: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html#linux-and-macos
- Windows: https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html#windows

Running the env installation will also take 5-10 mins, and it needs to download and extract pyrosetta. By default, it install the environment in the current directory. If you want to install it in a different directory, modify or remove the `prefix` in the: `micromamba env create` line.

### Quick start

PatchMAN can be run with:

```
micromamba activate ./patchman   # default the env is created inside the directory. 
python3 PatchMAN_protocol.py <arguments> RECEPTOR PEPTIDE
```

where RECEPTOR is a PDB file and PEPTIDE is the peptide sequence to be docked.
The peptide can contain post-translational modifications, denoted by Rosetta standards, e.g. `[SER:phosphorylated]`. The available PTMs can be listed with

```ls <env_path>/lib/python3.11/site-packages/pyrosetta/database/chemical/residue_type_sets/fa_standard/patches```

:exclamation: Note that the protocol script is set up to use Slurm job scheduler. Using an other type of scheduler needs editing of the `config.ini` file and the `.sh`  and fpd.py files in the `bin/` directory. Unfortunately, we cannot help with that.

#### Test run
A test run of PatchMAN can be performed on the 1ssh.pdb in the `test/` directory. Turning off receptor backbone minimization for testing purposes decreases runtime:

```
python3 PatchMAN_protocol.py 1ssh.pdb -w test/ EGPPPAMPARPT
```

### Running parameters
```
Focus/mask/hotspot:
-s mask/focus PDB file (type: pdb file, Default: None)
-l mask/focus residue list, instead of PDB file
-f focus mode: use the file specified as the focus (Default: False - mask mode)
-o hotspot mode, with a few residues in the binding site (Default: False)

# MASTER parameters
-c master cutoff (Default: 1.5)

# FPD related parameters
-m minimize receptor backbone (default: false)
-a native file, for benchmarking (type: pdb file, Default: None)
-t number of structures to generate (Default: 1)
-u cluster radius (Default: 2.0)

-w working directory (Default: current directory)
-p step to start from (Default: 1, 1: split to motifs, 2: prepack receptor, 3: run MASTER,
                                4: extract templates,  5: FlexPepDock, 6: clustering and finalizing,
                                example: 4-6)
                                
-v verbose mode, print more information about the job
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
