#!/bin/bash
#SBATCH --nice=8000
#SBATCH	--time=100:00:00
#SBATCH	--kill-on-invalid-dep=yes
#SBATCH	--mem-per-cpu=1600
#SBATCH	--get-user-env
#SBATCH --ntasks=150
#SBATCH --module="openmpi/2.1.6"

module load singularity openmpi/2.1.6
ls *0001.pdb > input_list
export HWLOC_HIDE_ERRORS=1 # otherwise MPI crashes on multiple nodes

n_templates=`wc -l input_list | gawk '{print $1}'`
min_nstruct=`gawk "BEGIN {print $n_templates * 10}"`
max_nstruct=20000

# Set default parameters
nstruct=''
native=''
min_rec_bb=false

usage() {
	cat <<-USAGE
	Usage: ${0##*/} [opts]
	This script controls the FlexPepDock refinement part of the protocol. The following arguments can be set:

	-t number of FlexPepDock decoys (nstruct) to generate. Default: 1
	-m use receptor backbone minimization? Approximately doubles the runtime. Default: false
	-u unboundrot. Add unbound rotamers for the receptor. Default: none
	-a native structure for comparison. Default: none

	USAGE
}

while getopts :hu:t:a:m: opt; do
	case $opt in
		h)
			usage
			exit 0
			;;
		u)
			unboundrot=$OPTARG
			;;
		t)
			nstruct=$OPTARG
			;;
		m)
			min_rec_bb=$OPTARG
			;;
		a)
			native=$OPTARG
			;;
		\?)
			echo "Invalid option: $OPTARG" >&2
			exit 1
			;;
		:)
			echo "Requires an argument $OPTARG" >&2
			exit 1
			;;
	esac
done
shift "$((OPTIND-1))"

if [ -z $nstruct ]; # if it not has been set manually
then
	echo "Nstruct needs to be set automatically"
	if [ $min_nstruct -lt $max_nstruct ];
	then
		nstruct=10;
	elif [ $n_templates -gt $max_nstruct ];
	then
		nstruct=1;
	else
		nstruct=`gawk "BEGIN {print $max_nstruct / $n_templates}" | cut -d '.' -f 1`
	fi;
fi
echo "Refinement will generate $nstruct decoys per input structure"

args="-in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary \
-in:file:l input_list -scorefile score.sc -out:file:silent_struct_type binary -overwrite \
-out:file:silent decoys.silent -lowres_preoptimize -flexPepDocking:pep_refine -flexPepDocking:flexpep_score_only \
-ex1 -ex2aro -use_input_sc -unboundrot "$unboundrot" -min_receptor_bb "$min_rec_bb" -nstruct $nstruct"

if [ -z "$native" ]
then
  mpirun ${ROSETTA_BIN}/FlexPepDocking.mpiserialization.linuxgccrelease $args -native "$native"
else
  mpirun ${ROSETTA_BIN}/FlexPepDocking.mpiserialization.linuxgccrelease $args
fi
