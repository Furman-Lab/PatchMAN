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
#           Created by Furman Lab at HUJI, 2023.          #
###########################################################


die() {
	echo >&2 -e "\nERROR: $@\n"
	exit 1
}

# if count_atom_lines greather than 0, then the file is a PDB file, return true
validate_pdb() {
	# if it is
	if [[ ! -r $1 ]]; then
		die "$1 is not a readable file"
	fi

  count_atom_lines=$(grep -Ec "^ATOM  [ 0-9]{5} [A-Z0-9 ']{4}[A-Z ][A-Z0-9 ]{3} [A-Z ][ 0-9]{4}[A-Z ] {4}[0-9. -]{8}[0-9. -]{8}[0-9. -]{7}[0-9 .]{6}[ 0-9.]{6}" $1)

  if [[ $count_atom_lines -gt 0 ]]; then
    return 0
  else
  	die "$1 is not a valid PDB file"
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

# copy the file to the working directory and return the absolute path of the file
prepare_pdb(){
		validate_pdb $1 # this dies if the file is not a valid PDB file
		cp $1 .
		new_file=$(readlink -f $(basename "$1"))

		echo $new_file
}

# print the jobid if in verbose mdoe
print_jobid(){
	[[ $verbose ]] || echo "DEBUG| $1 JOBID: " $2
}

# Defaults
export work_dir=$(pwd)
job_name="PatchMAN_JOB"
cluster_radius="2.0"
min_rec_bb="true"
nstruct=''
master_cutoff="1.5"
step_from=1
step_to=6
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
					-o hotspot mode, only the residues in focus will be used as patch centers (Default: false)
					-s PDB file with residues that should or should not be in the binding site (type: pdb file, Default: None)
					-l list of residues that should or should not be in the binding site (type: string, Default: None, format: 'A23,A12')
					-f focus mode - the residues should from -l or -s should form the binding site (type: boolean, Default: False, masking mode)
					-p steps to run between (Default: 1-6,
									1: split to motifs, 2: prepack receptor, 3: run MASTER,
									4: extract templates,  5: FlexPepDock, 6: clustering and finalizing)
					-b benchmark mode, use the benchmark mode of FlexPepDock with the file provided (Default: None)
	USAGE
}

while getopts hvw:g:c:t:fs:n:m:p:f:a:orb:l: opt; do
	case $opt in
		h)
			usage
			exit 0
			;;
		a)
			native=$(readlink -f $OPTARG)
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
			IFS="-" read step_from step_to <<< "$OPTARG"
			[ -z "$step_from" ] && step_from=1
			[ -z "$step_to" ] && step_to=6
			;;
		r)
			fpd_args="$fpd_args -r "
			;;
		s)
			mask_focus=$(readlink -f $OPTARG)
			;;
		f)
			focus_mask_args="$focus_mask_args -f"
			;;
		l)
			resi_list=$OPTARG
			focus_mask_args="$focus_mask_args -l $resi_list"
			;;
		o)
			hotspot_mode=1
			focus_mask_args="$focus_mask_args -o "
			;;
		b)
			benchmark=$OPTARG
			;;
		v)
			verbose=1
			args=" -v "
			;;
		n)
			job_name=$(echo $OPTARG | sed -r 's/[\t\n ]+/_/g')
			;;
		\?)
			echo "Invalid option: $OPTARG" >&2
			 1
			;;
		:)
			echo "Requires $OPTARG" >&2
			exit 1
			;;
	esac
done
shift "$((OPTIND-1))"

# Get absolute path of the input receptor
receptor=$(readlink -f $1)


###########################################################
# Set protocol root based on  of this script.
export PROTOCOL_ROOT=$(dirname $(realpath ${BASH_SOURCE}))
root_dir=$(echo "${work_dir}" | grep -oE '^/[^/]+')
source "${PROTOCOL_ROOT}/.env" ",${root_dir}"
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

# Create output directory if does not exist, and cd into it
mkdir -p $work_dir
pushd $work_dir > /dev/null

############# VALIDATE INPUT #############

# Prepare receptor
# If validate_pdb returns 0, then the receptor is valid
receptor=$(prepare_pdb $receptor)
receptor_base=$(basename "$receptor")

# Prepare input peptide
pep_sequence_to_validate=$(echo "$2" |  sed 's/\[[A-Z]{3,4}\:[a-z]+\]//g' -E) # remove PTMs for validation of the rest of the peptide
[[ "$pep_sequence_to_validate" =~ ^[ARNDCEQGHILKMFPSTWYV]+$ ]] || die "Not a peptide sequence: '$2'" # modified for PTM
pep_sequence="$2"

# If both a pdb file or a list of residues are provided, then raise an error
# Check if files are readable ones and if not, raise an error
# If they are okay, copy to the working directory, get their path and add to the focus_mask_args
if [[ (-n "$focus_list" && -n "$mask_focus") ]]; then
	die "Provide either PDB file or list for mask/focus residues!"
elif [[ -n $mask_focus ]]; then
	mask_focus=$(prepare_pdb $mask_focus)
	focus_mask_args="$focus_mask_args -s $mask_focus"
fi

if [[ -n $native ]]; then
	native=$(prepare_pdb $native)
	fpd_args="$fpd_args -a $native"
fi

if [[ -n $benchmark ]]; then
	# grep the first 4 letters of the input pdb from the benchmark file and write a new one, splitting the line by '|'
	grep -m 1 -Ei "${receptor_base:0:4}" $benchmark | sed 's/|/\n/g' > benchmark_file
	fpd_args="$fpd_args -b"
fi

############### PREPARE JOB ###############

# The protocol can only handle one chain. If more than one chains are in the receptor, throw an error.
chain_ids=$(grep '^ATOM' $receptor | cut -c 22 | sort | uniq | wc -l)

if [[ $chain_ids -lt 1 ]]; then 
	die "More than one chain is provided for the receptor. The protocol can only handle one chain. Rename your chains to run PatchMAN"
fi

# Set pdb filenames
rec_name=`echo ${receptor_base::-4}`
clean_rec="$rec_name.clean.pdb"
ppkrec=`echo ${receptor_base::-4}'.clean.ppk.pdb'`
[[ $verbose ]] || echo "DEBUG: " $clean_rec $rec_name $ppkrec

# Step 1: Split to motifs
if [[ 1 -ge $step_from && 1 -le $step_to ]]
then
	$PYTHON ${BIN_DIR}/split_to_motifs.py -i "$receptor" $args $focus_mask_args ||
	"Splitting receptor to patches was not successful, aborting"

	ls ???'_'$rec_name'.pdb' > motif_list
	$MASTER/createPDS --type query --pdbList motif_list >& /dev/null # remove the long stdout
	echo "MASTER pds files were created for all motifs"
else
	echo "Skipping step 1: Splitting receptor to motifs"
fi

# This will be used by several steps
n_searches=$(wc -l motif_list | gawk '{print $1}')

# Step 2: Prepack receptor
if [[ 2 -ge $step_from && 2 -le $step_to ]]
then
	prepack_receptor_jid=$(sbatch --job-name=prepack_receptor --get-user-env --time=90:00:00\
	                --mem=1600m ${BIN_DIR}/prepare_input.sh $clean_rec | awk '{print $NF}')
else
	echo "Skipping step 2: Prepacking receptor"
fi

# Step 3: Run MASTER
if [[ 3 -ge $step_from && 3 -le $step_to ]]
then
	prepack_receptor_jid=$(zero_jobid $prepack_receptor_jid)
	print_jobid "PREPACK" $prepack_receptor_jid
	run_master_jid=$(sbatch --dependency=afterok:"${prepack_receptor_jid}" --array=0-"$n_searches"%50 ${BIN_DIR}/run_master.sh $master_cutoff | awk '{print $NF}')
else
	echo "Skipping step 3: Running MASTER"
fi

# Step 4: Extract templates
if [[ 4 -ge $step_from && 4 -le $step_to ]]
then
	run_master_jid=$(zero_jobid $run_master_jid)
	print_jobid "MASTER" $run_master_jid
	extract_templates_jid=$(sbatch --array=0-"$n_searches"%50 --dependency=afterany:"${run_master_jid}" ${BIN_DIR}/run_extract_templates.sh \
	                    "$pep_sequence" "$ppkrec" "$focus_mask_args" | awk '{print $NF}')
else
	echo "Skipping step 4: Extracting templates"
fi

# Step 5: FPD
if [[ 5 -ge $step_from && 5 -le $step_to ]]
then
	extract_templates_jid=$(zero_jobid $extract_templates_jid)
	print_jobid "EXTRACT TEMPLATES" $extract_templates_jid

	if [ ! -z "$nstruct" ]
	then
		fpd_args="$fpd_args -t $nstruct"
	fi
	echo Running FPD
	fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd $ADD_SBATCH \
					fpd.sh -u "$clean_rec" -m "$min_rec_bb" $fpd_args | awk '{print $NF}')
else
	echo "Skipping step 5: FlexPepDock"
fi


# Step 6: Clustering & Step 6: Finalizing
if [[ 6 -ge $step_from && 6 -le $step_to ]]
then
	print_jobid "FPD" $fpd_jid

	echo "Running steps: clustering and finalizing"
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
else
	echo "Skipping step 6: Clustering and finalizing"
fi




[ "$verbose" ] && {

	cat <<-JOBINFO
	------------------------------------------------
	$(date)
	Verbosity: $verbose
	Receptor is: $receptor
	Peptide sequence is: $pep_sequence
	Cluster radius: $cluster_radius
	Slurm job IDs: $run_master_jid $extract_templates_jid $fpd_jid $clustering_jid $finalize_jid
	Current working directory: $(pwd)
	------------------------------------------------
	JOBINFO

}

echo "SLURM_JID $finalize_jid"
