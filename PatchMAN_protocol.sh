#!/bin/bash

###########################################################
#                  The PatchMAN protocol                  #
#                                                         #
# The protocol splits the surface into patches, then      #
# searches PDB30 with these to find peptide fragments     #
# that can complement them. These are extracted and used  #
# as templates for docking the peptide.                   #
# This script sends jobs asynchronously via Slurm, using  #
# Singularity containers.                                 #
#                                                         #
#           Created by Furman Lab at HUJI, 2022.          #
###########################################################


die() {
	echo >&2 -e "\nERROR: $@\n"
	 1
}

# if count_atom_lines greather than 0, then the file is a PDB file, return true
validate_pdb() {
  count_atom_lines=$(grep -Ec "^ATOM  [ 0-9]{5} [A-Z0-9 ']{4}[A-Z ][A-Z0-9 ]{3} [A-Z ][ 0-9]{4}[A-Z ] {4}[0-9. -]{8}[0-9. -]{8}[0-9. -]{8}[0-9 .]{6}[ 0-9.]{6} {5}[A-Z ]{2}.{0,6}" $1)
  if [[ $count_atom_lines -gt 0 ]]; then
    return 0
  else
    return 1
  fi
}

# if jobid is empty, set it to 0
zero_jobid() {
	if [[ -z "$1" ]]; then
		echo -1
	else
		echo $1
	fi
}

[ -d $PROTOCOL_ROOT ] || die "Protocol root directory is not a directory: $PROTOCOL_ROOT"

###########################################################
# Set protocol root based on  of this script.
export PROTOCOL_ROOT=$(dirname $(realpath ${BASH_SOURCE}))
source ${PROTOCOL_ROOT}/.env $PROTOCOL_ROOT
export BIN_DIR=${PROTOCOL_ROOT}/bin
export PATH=.:${BIN_DIR}:${PATH}
export PYTHONPATH="" # messes up python packages inside the  otherwise
[ -d $PROTOCOL_ROOT ] || die "Protocol root directory is not a directory: ${PROTOCOL_ROOT}"
###########################################################
if [[ "$VIRTUAL_ENV" == '' ]]
then
        export PYTHON=$(echo $PYTHON | sed "s#PROTOCOL_ROOT#${PROTOCOL_ROOT}#g")
else
        . $VIRTUAL_ENV/bin/activate || die "No virtual  detected. Please install it  by: virtualenv .venv && . .venv/bin/activate && pip install -r requirements.txt"
        export PYTHON="python3 "
fi
###########################################################
#module try-load openmpi/2.1.6

# Defaults
work_dir=$(pwd)
job_name="PatchMAN_JOB"
cluster_radius="2.0"
min_rec_bb="true"
nstruct=''
master_cutoff="1.5"
step=1
verbose=False

usage() {
	cat <<-USAGE
	Usage: ${0##*/} [opts] RECEPTOR PEPTIDE_SEQUENCE
	PatchMAN performs search on existing  monomers and complexes with structural motifs extracted from the query receptor and extract complementary fragments to be used as templates for peptide-protein interactions.
		RECEPTOR: PDB file with the receptor 
		PEPTIDE_SEQUENCE can include modified residues in "GFK[SER:phosphorylated]RAD" format.
        	-m minimize receptor backbone (default: false)
					-g log file (Default is stdout)
					-e error log file (Default is stderr)
					-n job name (Default: PatchMAN_JOB)
	        -v print information about the job
	        -w working directory (Default: current directory)
	        -t number of structures to generate (Default: 1)
	        -c master cutoff (Default: 1.5)
	        -s mask file with resides not in the binding site (type: pdb file, Default: None)
	        -f focus mask, with residues that are in the binding site (type: pdb file, Default: None)
	        -p step to start from (Default: 1, 1: split to motifs, 2: prepack receptor, 3: run MASTER,
	        												4: extract templates,  5: FlexPepDock, 6: clustering and finalizing)
	USAGE
}

while getopts :hvw:g:c:t:f:s:n:m:p:f:a: opt; do
	case $opt in
		h)
			usage
			exit 0
			;;
		a)
			native=$OPTARG
			;;
		g)
			logs_dir=$OPTARG
			;;
		c)
		  master_cutoff=$OPTARG
		  ;;
		w)
			export work_dir=$(realpath $OPTARG)
			;;
		t)
			nstruct=$OPTARG
			;;
		m)
			min_rec_bb=$OPTARG
			;;
		p)
			step=$OPTARG
			;;
		s)
			mask=$(readlink -f $OPTARG)
			;;
		f)
			focus=$(readlink -f $OPTARG)
			;;
		v)
			verbose=True
			;;
		n)
			job_name=$(echo $OPTARG | sed -r 's/[\t\n ]+/_/g')
			;;
		\?)
			echo "Invalid option: $OPTARG" >&2
			 1
			;;
		:)
			echo "Requires  $OPTARG" >&2
			exit 1
			;;
	esac
done
shift "$((OPTIND-1))"

[ -r "$1" ] || die "Receptor is not a readable file: $1"
pep_sequence_to_validate=$(echo "$2" |  sed 's/\[[A-Z]{3,4}\:[a-z]+\]//g' -E) # remove PTMs for validation of the rest of the peptide
[[ "$pep_sequence_to_validate" =~ ^[ARNDCEQGHILKMFPSTWYV]+$ ]] || die "Not a peptide sequence: $2" # modified for PTM

############### PREPARE JOB ###############
# Creating a directory for the job and copying inputs to it
receptor=$(readlink -f $1)
pep_sequence="$2"

# Create output directory if does not exist, and cd into it
mkdir -p $work_dir
pushd $work_dir > /dev/null || 

# Prepare receptor
# if validate_pdb returns 0, then the receptor is valid
validate_pdb $receptor || die "Receptor is not a valid PDB file: $receptor"
cp $receptor .
receptor=$(readlink -f $(basename "$receptor"))
receptor_base=$(basename "$receptor")

# Rename all receptor chains to one, otherwise the protocol crashes
chain_id=$(grep -m 1 '^ATOM' $receptor | cut -c 22)
sed -Ei "s/^(ATOM.{17})[A-Z]/\1${chain_id}/" $receptor
sed '/^TER/d' -i $receptor

# Prepare (focus) mask if provided
if  [[ -f "$mask" ]]
then
  # if validate_pdb returns 0, then the mask is valid
  validate_pdb $mask || die "Mask is not a valid PDB file: $mask"
  cp $mask .
  mask=$(readlink -f $(basename "$mask"))
elif [[ -f "$focus" ]]
then
	# if validate_pdb returns 0, then the mask is valid
	validate_pdb $focus || die "Mask is not a valid PDB file: $focus"
	cp $focus .
	focus=$(readlink -f $(basename "$focus"))
fi

# Set pdb filenames
rec_name=`echo ${receptor_base::-4}`
clean_rec="$rec_name.clean.pdb"
ppkrec=`echo ${receptor_base::-4}'.clean.ppk.pdb'`
echo "DEBUG| " $clean_rec $rec_name $ppkrec

# Step 1: Split to motifs
if [[ $step -le 1 ]]
then
	# if mask or focus is provided, append to args
	if [[ -f "$mask" ]]
	then
		args="-m $mask"
	elif [[ -f $focus ]]; then
		args="-f $focus"
	else
		args=""
	fi

	# if verbose, append to args
	if [[ $verbose == True ]]
	then
		args="$args -v"
	fi

	$PYTHON ${BIN_DIR}/split_to_motifs.py -i "$receptor" $args

	ls ???'_'$rec_name'.pdb' > motif_list
	$MASTER/createPDS --type query --pdbList motif_list >& /dev/null # remove the long stdout
	echo "MASTER pds files were created for all motifs"
	n_searches=$(wc -l motif_list | gawk '{print $1}')
fi

# Step 2: Prepack receptor
if [[ $step -le 2 ]]
then
	prepack_receptor_jid=$(sbatch --job-name=prepack_receptor --get-user-env --time=90:00:00\
	                --mem=1600m ${BIN_DIR}/prepare_input.sh $clean_rec | awk '{print $NF}')
fi

# Step 3: Run MASTER
if [[ $step -le 3 ]]
then
	prepack_receptor_jid=$(zero_jobid $prepack_receptor_jid)
	if [[ $verbose == True ]]
	then
		echo "DEBUG| PREPACK JOBID: " $prepack_receptor_jid
	fi
	run_master_jid=$(sbatch --dependency=afterok:"${prepack_receptor_jid}" --array=0-"$n_searches"%50 ${BIN_DIR}/run_master.sh $master_cutoff | awk '{print $NF}')
fi

# Step 4: Extract templates
if [[ $step -le 4 ]]
then
	run_master_jid=$(zero_jobid $run_master_jid)
	if [[ $verbose == True ]]
	then
		echo "DEBUG| MASTER JOBID: " $run_master_jid
	fi
	extract_templates_jid=$(sbatch --array=0-"$n_searches"%50 --dependency=afterany:"${run_master_jid}" ${BIN_DIR}/run_extract_templates.sh \
	                    "$pep_sequence" "$ppkrec" | awk '{print $NF}')
fi

# Step 5: FPD
if [[ $step -le 5 ]]
then
	extract_templates_jid=$(zero_jobid $extract_templates_jid)
	if [[ $verbose == True ]]
	then
		echo "DEBUG| EXTRACT TEMPLATES JOBID: " $extract_templates_jid
	fi

	extra_args=''
	if [ -z "$native" ]
	then
		extra_args="-a $native"
	fi

	if [ -z "$nstruct" ]
	then
		extra_args="$extra_args -t $nstruct"
	fi

	fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd \
					fpd.sh -u "$clean_rec" -m "$min_rec_bb" $extra_args | awk '{print $NF}')
fi

# Step 6: Clustering & Step 6: Finalizing
if [[ $step -le 6 ]]
then
	if [[ $verbose == True ]]
	then
		echo "DEBUG| FPD JOBID: " $fpd_jid
	fi

	echo "Running final steps"
	fpd_jid=$(zero_jobid $fpd_jid)
	clustering_jid=$(sbatch \
					--job-name=clustering \
					--nice=8000 \
			--chdir=$(pwd) \
			--dependency=aftercorr:${fpd_jid} \
			--kill-on-invalid-dep=yes \
			--get-user-env \
			cluster.sh $pep_sequence $cluster_radius | awk '{print $NF}')

	clustering_jid=$(zero_jobid $clustering_jid)
	finalize_jid=$(sbatch \
					--job-name=finalize \
					--nice=8000 \
			--chdir=$(pwd) \
			--dependency=aftercorr:${clustering_jid} \
			--kill-on-invalid-dep=yes \
			--mem-per-cpu=2000 \
			--get-user-env \
			finalize.sh | awk '{print $NF}')
fi






[ "$verbose" ] && {

	cat <<-JOBINFO
	------------------------------------------------
	$(date)
	Receptor is: $receptor
	Peptide sequence is: $pep_sequence
	Cluster radius: $cluster_radius
	Slurm job IDs: $run_master_jid $extract_templates_jid $fpd_jid $clustering_jid $finalize_jid
	Current working directory: $(pwd)
#	Notification script: $([ "$notify_script" ] && echo $notify_script)
	------------------------------------------------
	JOBINFO

}

echo "SLURM_JID $finalize_jid"
