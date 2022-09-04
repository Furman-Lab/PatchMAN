#!/bin/bash
#SBATCH --nice=8000
#SBATCH	--time=100:00:00
#SBATCH	--kill-on-invalid-dep=yes
#SBATCH	--mem-per-cpu=1600
#SBATCH	--get-user-env
#SBATCH --ntasks=150
#SBATCH --module="openmpi/2.1.6"

#module load openmpi/2.1.6
ls *0001.pdb > input_list;
mpirun $ROSETTA_BIN/FlexPepDocking.mpiserialization.linuxgccrelease -in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary \
-out:file:silent decoys.silent -lowres_preoptimize -flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only \
-ex1 -ex2aro -use_input_sc -unboundrot "$1" -min_receptor_bb "$2" 
