#!/bin/bash
#SBATCH --nice=8000
#SBATCH	--time=100:00:00
#SBATCH	--kill-on-invalid-dep=yes
#SBATCH	--mem-per-cpu=1600
#SBATCH	--get-user-env
#SBATCH --ntasks=150

ls *0001.pdb > input_list;
n_templates=`wc -l input_list | gawk '{print $1}'`;
min_nstruct=`gawk "BEGIN {print $n_templates * 10}"`;
max_nstruct=20000;
if [ $min_nstruct -lt $max_nstruct ];
then
  nstruct=10;
elif [ $n_templates -gt $max_nstruct ];
then
  nstruct=1;
else
  nstruct=`gawk "BEGIN {print $max_nstruct / $n_templates}" | cut -d '.' -f 1`;
fi;
echo $nstruct;

if (($# == 3))
then
  mpirun $ROSETTA_BIN/FlexPepDocking.mpiserialization.linuxgccrelease -in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary \
  -out:file:silent decoys.silent -lowres_preoptimize -flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only \
  -ex1 -ex2aro -use_input_sc -unboundrot "$1" -min_receptor_bb "$2" -native "$3" -nstruct "$nstruct"
elif (($# == 2))
then
  mpirun $ROSETTA_BIN/FlexPepDocking.mpiserialization.linuxgccrelease -in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary \
  -out:file:silent decoys.silent -lowres_preoptimize -flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only \
  -ex1 -ex2aro -use_input_sc -unboundrot "$1" -min_receptor_bb "$2" -nstruct "$nstruct"
else
  echo "Missing arguments for refinement"
  echo "The provided arguments are:"
  echo $@
fi
