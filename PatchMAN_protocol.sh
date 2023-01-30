#!/bin/bash
# The PatchMAN protocol.
# It sends jobs asynchronously via Slurm. It is very important that this scripts
# returns as fast as possible, as it is called periodically by the daemon, so all actions here needs to
# to be fast - no copying large files, running external scripts synchronously etc.

die() {
	echo >&2 -e "\nERROR: $@\n"
	exit 1
}

# if count_atom_lines greather than 0, then the file is a PDB file, return true
validate_pdb() {
  count_atom_lines=$(grep -Ec "^ATOM  [ 0-9]{5} [A-Z0-9 ']{4}[A-Z ][A-Z0-9 ]{3} [A-Z ][ 0-9]{4}[A-Z ] {4}[0-9. -]{8}[0-9. -]{8}[0-9. -]{8}[0-9 .]{6}[ 0-9.]{6} {9}[A-Z ]{2}[A-Z ]{0,2}" $1)
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

export PROTOCOL_ROOT=$(dirname $(realpath $BASH_SOURCE)) # changed for development purposes but could probably stay like that
export BIN_DIR=${PROTOCOL_ROOT}/bin

export ROSETTA_DB=/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/main/database
export ROSETTA_BIN=/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/main/source/bin
export ROSETTA_TOOLS=/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/tools
export MASTER=/vol/ek/share/master_forPatchMAN

export PATH=.:${BIN_DIR}:/vol/ek/share/labscripts:/vol/ek/share/bin:usr/local/bin:/usr/bin:/bin:/usr/X11R6/bin:/usr/lib/mh:/etc/alternatives/slurm/bin:/vol/slurm/picasso/bindir/bin/
export VIRTUAL_ENV=/cs/labs/fora/projects/autopeptidb/staging/venv_PatchmanProtocol
#export VIRTUAL_ENV=/cs/labs/fora/projects/autopeptidb/staging/virtualenv-linux


#activating virtual env (needed for various python libraries used in the protocol)
. $VIRTUAL_ENV/bin/activate || die "No virtual environment detected. Please install it first by: virtualenv .venv && . .venv/bin/activate && pip install -r requirements.txt"
export PYTHONPATH=''
module load openmpi/2.1.6

# Defaults
work_dir=$(pwd)
job_name="PatchMAN_JOB"
cluster_radius="2.0"
min_rec_bb="true"
nstruct=1
master_cutoff="1.5"
step=1

usage() {
	cat <<-USAGE
	Usage: ${0##*/} [opts] RECEPTOR PEPTIDE_SEQUENCE
	TODO: describe

		TODO: replace this with real options
        	-m minimize receptor backbone (default: false)
					-g log file (Default is stdout)
					-e error log file (Default is stderr)
					-n job name (Default: PatchMAN_JOB)
	        -v print information about the job
	        -w working directory (Default: current directory)
	        -t number of structures to generate (Default: 1)
	        -c master cutoff (Default: 1.5)
	        -s mask file (Default: None)
	        -p step to start from (Default: 1, 1: split to motifs, 2: prepack receptor, 3: run MASTER,
	        												4: extract templates,  5: FlexPepDock, 6: clustering and finalizing)


	USAGE
}

while getopts :hvw:g:c:t:f:s:n:m:p: opt; do
	case $opt in
		h)
			usage
			exit 0
			;;
		g)
			logs_dir=$OPTARG
			;;
		c)
		  master_cutoff=$OPTARG
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
		p)
			step=$OPTARG
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
[[ "$pep_sequence_to_validate" =~ ^[ARNDCEQGHILKMFPSTWYV]+$ ]] || die "Not a peptide sequence: $2" # modified for PTM

# Creating a directory for the job and copying inputs to it
receptor=$(readlink -f $1)
pep_sequence=$2

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
if  [[ -f "$mask" ]]
then
  # if validate_pdb returns 0, then the mask is valid
  validate_pdb $mask || die "Mask is not a valid PDB file: $mask"
  cp $mask .
  mask=$(readlink -f $(basename "$mask"))
fi

# Set pdb filenames
rec_name=`echo ${receptor_base::-4}`
clean_rec="$rec_name.clean.pdb"
ppkrec=`echo ${receptor_base::-4}'.clean.ppk.pdb'`
echo "DEBUG| " $clean_rec $rec_name $ppkrec

# Step 1: Split to motifs
if [[ $step -le 1 ]]
then
	python ${BIN_DIR}/split_to_motifs.py "$receptor" "$mask"
	ls ???'_'$rec_name'.pdb' > motif_list
	$MASTER/createPDS --type query --pdbList motif_list
	n_searches=$(wc -l motif_list | gawk '{print $1}')
fi

# Step 2: Prepack receptor
if [[ $step -le 2 ]]
then
	prepack_receptor_jid=$(sbatch --job-name=prepack_receptor --get-user-env --time=90:00:00\
	                --mem=1600m prepack_receptor.sh $clean_rec | awk '{print $NF}')
fi


# Step 3: Run MASTER
if [[ $step -le 3 ]]
then
	prepack_receptor_jid=$(zero_jobid $prepack_receptor_jid)
	run_master_jid=$(sbatch --dependency=afterok:"${prepack_receptor_jid}" --array=0-"$n_searches"%50 run_master.sh $master_cutoff | awk '{print $NF}')
fi

# Step 4: Extract templates
if [[ $step -le 4 ]]
then
	run_master_jid=$(zero_jobid $run_master_jid)
	extract_templates_jid=$(sbatch --array=0-"$n_searches"%50 --dependency=after:"${run_master_jid}" run_extract_templates.sh \
	                    "$pep_sequence" "$ppkrec" | awk '{print $NF}')
fi

# Step 5: FPD
if [[ $step -le 5 ]]
then
	echo 'DEBUG|FPD'
	extract_templates_jid=$(zero_jobid $extract_templates_jid)
	fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd \
						fpd.sh "$clean_rec" "$min_rec_bb" "$nstruct" | awk '{print $NF}')
fi

# Step 6: Clustering & Step 6: Finalizing
if [[ $step -le 6 ]]
then
	echo 'DEBUG|submitting clustering and finalize'
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
