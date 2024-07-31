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

  count_atom_lines=$(grep -Ec "^ATOM  [ 0-9]{5} [A-Z0-9 ']{4}[A-Z ][A-Z0-9 ]{3} [A-Z ][ 0-9]{4}[A-Z ] {4}[0-9. -]{8}[0-9. -]{8}[0-9. -]{8}[0-9 .]{6}[ 0-9.]{6}[ A-Z0-1]{9}[0-9A-Z -]{2,4}[A-Z ]{0,2}" $1)
  if [[ $count_atom_lines -gt 0 ]]; then
    return 0
  else
    return 1
  fi
}

[ -d $PROTOCOL_ROOT ] || die "Protocol root directory is not a directory: $PROTOCOL_ROOT"

export PROTOCOL_ROOT=$(dirname $(realpath $BASH_SOURCE)) # changed for development purposes but could probably stay like that
export BIN_DIR=${PROTOCOL_ROOT}/bin

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
cluster_radius="2.0"
min_rec_bb="false"
master_cutoff="1.5"
step_from=1
step_to=6
verbose=False

usage() {
	cat <<-USAGE
	Usage: ${0##*/} [opts] RECEPTOR PEPTIDE_SEQUENCE

	PatchMAN performs search on existing  monomers and complexes with structural motifs extracted from the query receptor and extract complementary fragments to be used as templates for peptide-protein interactions.

		RECEPTOR: PDB file with the receptor
		PEPTIDE_SEQUENCE

					-m minimize receptor backbone (default: false)
					-g log file (Default is stdout)
					-e error log file (Default is stderr)
					-n job name (Default: PatchMAN_JOB)
					-v print information about the job (debug)
					-w working directory (Default: current directory)
					-c master cutoff (Default: 1.5)
					-a native structure (PDB format) for comparison
				        -p steps to run between (Default: 1-6,
									1: split to motifs, 2: prepack receptor, 3: run MASTER,
									4: extract templates,  5: FlexPepDock, 6: clustering and finalizing)
	USAGE
}

while getopts hvw:c:m:a:p: opt; do
	case $opt in
		h)
			usage
			exit 0
			;;
		a)
			native=$(readlink -f $OPTARG)
			;;
		c)
			master_cutoff=$OPTARG
			;;
		w)
			export work_dir=$(realpath $OPTARG)
			;;
		p)
			IFS="-" read step_from step_to <<< "$OPTARG"
			[ -z "$step_from" ] && step_from=1
			[ -z "$step_to" ] && step_to=6
			;;
		m)
			min_rec_bb=$OPTARG
			;;
		v)
			verbose=1
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
receptor=$(prepare_pdb "$receptor")
receptor_base=$(basename "$receptor")

# Prepare input peptide
pep_sequence="$2"
[[ "$pep_sequence" =~ ^[ARNDCEQGHILKMFPSTWYV]+$ ]] || die "Not a peptide sequence: '$2'"

############### PREPARE JOB ###############

# The protocol can only handle one chain. If more than one chains are in the receptor, throw an error.
chain_ids=$(grep '^ATOM' "$receptor" | cut -c 22 | sort | uniq | wc -l)

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
	$PYTHON ${BIN_DIR}/split_to_motifs.py "$receptor" ||
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
	prepack_receptor_jid=$(sbatch --job-name=prep_input --get-user-env --time=90:00:00\
                --mem=1600m prepare_input.sh $clean_rec | awk '{print $NF}')
else
	echo "Skipping step 2: Prepacking receptor"
fi

# Step 3: Run MASTER
if [[ 3 -ge $step_from && 3 -le $step_to ]]
then
	prepack_receptor_jid=$(zero_jobid $prepack_receptor_jid)
	print_jobid "PREPACK" $prepack_receptor_jid
	run_master_jid=$(sbatch --dependency=afterany:"${prepack_receptor_jid}" --array=0-"$n_searches"%50 run_master.sh "$master_cutoff" | awk '{print $NF}')
else
	echo "Skipping step 3: Running MASTER"
fi

# Step 4: Extract templates
if [[ 4 -ge $step_from && 4 -le $step_to ]]
then
	run_master_jid=$(zero_jobid $run_master_jid)
	print_jobid "MASTER" "$run_master_jid"
	extract_templates_jid=$(sbatch --array=0-"$n_searches"%50 --dependency=afterok:"${run_master_jid}" run_extract_templates.sh \
                    "$pep_sequence" "$ppkrec" | awk '{print $NF}')
else
	echo "Skipping step 4: Extracting templates"
fi

# Step 5: FPD
if [[ 5 -ge $step_from && 5 -le $step_to ]]; then
	extract_templates_jid=$(zero_jobid $extract_templates_jid)
	if [[ ! -z "$native" ]]; then
       		fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd \
                 	fpd.sh "$clean_rec" "$min_rec_bb" "$native" | awk '{print $NF}')
	else
		fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd \
                  	fpd.sh "$clean_rec" "$min_rec_bb" | awk '{print $NF}')
	fi
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
