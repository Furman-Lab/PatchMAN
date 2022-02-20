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
	exit 1
}

# if count_atom_lines greather than 0, then the file is a PDB file, return true
validate_pdb() {
  count_atom_lines=$(grep -Ec "^ATOM  [ 0-9]{5} [A-Z0-9 ']{4}[A-Z ][A-Z0-9 ]{3} [A-Z ][ 0-9]{4}[A-Z ] {4}[0-9. -]{8}[0-9. -]{8}[0-9. -]{8}[0-9 .]{6}[ 0-9.]{6} {9}[A-Z ]{2}[A-Z ]{0,2}" $1)
  echo count $count_atom_lines
  if [[ $count_atom_lines -gt 0 ]]; then
    return 0
  else
    return 1
  fi
}

###########################################################
# Set protocol root based on the path of this script.
export PROTOCOL_ROOT=$(dirname $(realpath ${BASH_SOURCE}))
source ${PROTOCOL_ROOT}/.env $PROTOCOL_ROOT
export BIN_DIR=${PROTOCOL_ROOT}/bin
export PATH=.:${BIN_DIR}:${PATH}
export PYTHONPATH="" # messes up python packages inside the container otherwise
[ -d $PROTOCOL_ROOT ] || die "Protocol root directory is not a directory: ${PROTOCOL_ROOT}"
###########################################################
if [[ "$VIRTUAL_ENV" == '' ]]
then
	export PYTHON=$(echo $PYTHON | sed "s#PROTOCOL_ROOT#${PROTOCOL_ROOT}#g")
else
	. $VIRTUAL_ENV/bin/activate || die "No virtual environment detected. Please install it first by: virtualenv .venv && . .venv/bin/activate && pip install -r requirements.txt"
	export PYTHON="python3 "
fi

if [[ "$DB_PATH" == '' ]]
then
	export DB_PATH="${PROTOCOL_ROOT}/databases/master_clean/"
fi
###########################################################

# Defaults
work_dir=$(pwd)
job_name="PatchMAN_JOB"
cluster_radius="2.0"
min_rec_bb="true"
nstruct=1

usage() {
	cat <<-USAGE
	Usage: ${0##*/} [opts] RECEPTOR PEPTIDE_SEQUENCE
	PatchMAN performs search on existing protein monomers and complexes with structural motifs extracted from the query receptor and extract complementary fragments to be used as templates for peptide-protein interactions.

		RECEPTOR: PDB file with the receptor protein
		PEPTIDE_SEQUENCE can include modified residues in "GFK[SER:phosphorylated]RAD" format.

		-m minimize receptor backbone, longer runtime (true or false; default: true)
		-t number of refinement runs for FlexPepDock (default: 1)
		-s masked residues as a PDB structure (experimental, default: None)

		-w working directory (default: current directory)
		-o output directory for the ligands (default: working directory)
		-g log file (Default is	stdout)
		-e error log file (Default is stderr)
		-n job name
		-v verbose (default: false)
	USAGE
}


while getopts :hvw:g:t:f:j:s:n:m: opt; do
	case $opt in
		h)
			usage
			exit 0
			;;
		g)
			logs_dir=$OPTARG
			;;
		w)
			work_dir=$OPTARG
			;;
		t)
			nstruct=$OPTARG
			;;
		m)
			min_rec_bb=$OPTARG
			;;
		s)
			mask=$(readlink -f $OPTARG)
			;;
		v)
			verbose=True
			;;
		n)
			job_name=$(echo $OPTARG | sed -r 's/[\t\n ]+/_/g')
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

[ -r "$1" ] || die "Receptor is not a readable file: $1"
pep_sequence_to_validate=$(echo "$2" |  sed 's/\[[A-Z]{3,4}\:[a-z]+\]//g' -E) # remove PTMs for validation of the rest of the peptide
[[ "${pep_sequence_to_validate}" =~ ^[ARNDCEQGHILKMFPSTWYV]+$ ]] || die "Not a peptide sequence: $2" # modified for PTM

############### PREPARE JOB ###############
# Creating a directory for the job and copying inputs to it
receptor=$(readlink -f $1)
pep_sequence="$2"

pushd $work_dir > /dev/null

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

# Prepare mask if provided
if  [[ ! -n "$mask" ]]
then
  # if validate_pdb returns 0, then the mask file is an invalid PDB file
  validate_pdb $mask || die "Mask is not a valid PDB file: $mask"
  cp $mask .
  mask=$(readlink -f $(basename "$mask"))
fi


################ RUN JOB ################
# Step 1: Split to motifs
${PYTHON} ${PROTOCOL_ROOT}/bin/split_to_motifs.py "$receptor" "$mask"
clean_rec=`echo ${receptor_base::-4}`'.clean.pdb'
rec_name=`echo ${receptor_base::-4}`
ppkrec=`echo ${receptor_base::-4}'.clean.ppk.pdb'`
echo "DEBUG| " $clean_rec $rec_name $ppkrec
ls ???'_'$rec_name'.pdb' > motif_list
${MASTER}/createPDS --type query --pdbList motif_list

n_searches=$(wc -l motif_list | gawk '{print $1}')

# Step 2: Run MASTER
run_master_jid=$(sbatch --array=0-"$n_searches"%50 run_master.sh | awk '{print $NF}')


# Step2.5: Prepack receptor
prep_input_jid=$(sbatch --dependency=afterany:"${run_master_jid}" --job-name=prep_input --get-user-env --time=90:00:00\
                --mem=1600m prepare_input.sh $clean_rec | awk '{print $NF}')


# Step3: Extract templates
extract_templates_jid=$(sbatch --array=0-"$n_searches"%50 --dependency=afterok:"${prep_input_jid}" run_extract_templates.sh \
                    "$pep_sequence" "$ppkrec" | awk '{print $NF}')


# Step 4: FPD
fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd \
          fpd.sh "$clean_rec" "$min_rec_bb" "$nstruct" | awk '{print $NF}')


# Step 5: Clustering
clustering_jid=$(sbatch \
        --job-name=clustering \
        --nice=8000 \
		--chdir=$(pwd) \
		--dependency=aftercorr:${fpd_jid} \
		--kill-on-invalid-dep=yes \
		--get-user-env \
    cluster.sh $pep_sequence $cluster_radius | awk '{print $NF}')


# Step 6: Finalizing
finalize_jid=$(sbatch \
        --job-name=finalize \
        --nice=8000 \
		--chdir=$(pwd) \
		--dependency=aftercorr:${clustering_jid} \
		--kill-on-invalid-dep=yes \
		--mem-per-cpu=2000 \
		--get-user-env \
    finalize.sh | awk '{print $NF}')


[ "$verbose" ] && {

	cat <<-JOBINFO
	------------------------------------------------
	$(date)

	Receptor is: $receptor
	Peptide sequence is: $pep_sequence

	Cluster radius: $cluster_radius

	Slurm job IDs: $run_master_jid $extract_templates_jid $fpd_jid $clustering_jid $finalize_jid
	Current working directory: $(pwd)
	------------------------------------------------

	JOBINFO

}

echo "SLURM_JID $finalize_jid"
