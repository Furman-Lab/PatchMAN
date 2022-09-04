#!/bin/bash

die() {
	echo >&2 -e "\nERROR: $@\n"
	exit 1
}

[ -d $PROTOCOL_ROOT ] || die "Protocol root directory is not a directory: $PROTOCOL_ROOT"

export PROTOCOL_ROOT=$(dirname $(realpath $BASH_SOURCE)) # changed for development purposes but could probably stay like that
export BIN_DIR=${PROTOCOL_ROOT}/bin

export ROSETTA_DB=/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/main/database
export ROSETTA_BIN=/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/main/source/bin
export ROSETTA_TOOLS=/vol/ek/share/rosetta/rosetta_src_2019.14.60699_bundle/tools
export MASTER=/vol/ek/share/master_forPatchMAN

export PATH=.:${BIN_DIR}:/vol/ek/share/labscripts:/vol/ek/share/bin:${PATH}

export VIRTUAL_ENV=/sci/labs/fora/projects/autopeptidb/production/venv_PatchmanProtocol

#activating virtual env (needed for various python libraries used in the protocol)
. $VIRTUAL_ENV/bin/activate || die "No virtual environment detected. Please install it first by: virtualenv .venv && . .venv/bin/activate && pip install -r requirements.txt"


module load openmpi/2.1.6

# Defaults
work_dir=$(pwd)
job_name="PatchMAN_JOB"
cluster_radius="2.0"
min_rec_bb="false"
master_cutoff="1.5"


usage() {
	cat <<-USAGE
	Usage: ${0##*/} [opts] RECEPTOR PEPTIDE_SEQUENCE
	TODO: describe

		TODO: replace this with real options
		-o output directory for the ligands (Default: working directory)
		-p prefix for each ligand file extracted (Default: none)
        -s secondary structure (default: None)
        -m minimize receptor backbone (default: false)
		-g log file (Default is stdout)
		-e error log file (Default is stderr)
        -v print information about the job

	USAGE
}


while getopts :hvw:g:t:c:f:j:s:n:m: opt; do
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
			n_top_rots=$OPTARG
			;;
	  c)
      master_cutoff=$OPTARG
      ;;
    m)
			min_rec_bb=$OPTARG
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
[[ "$2" =~ ^[ARNDCEQGHILKMFPSTWYV]+$ ]] || die "Not a peptide sequence: $2"

# Creating a directory for the job and copying inputs to it
receptor=$(readlink -f $1)
pep_sequence=$2

# Rename all receptor chains to one, otherwise the protocol crashes
chain_id=$(grep -m 1 '^ATOM' $receptor | cut -c 22)
sed -Ei "s/^(ATOM.{17})[A-Z]/\1${chain_id}/" $receptor
sed -i "s/TER//g" $receptor
sed -i '/^$/d' $receptor

pushd $work_dir > /dev/null

cp $receptor .

receptor=$(readlink -f $(basename $receptor))
receptor_base=$(basename $receptor)

# Step 1: Split to motifs
python ${BIN_DIR}/split_to_motifs.py "$receptor"

clean_rec=`echo ${receptor_base::-4}`'.clean.pdb'
rec_name=`echo ${receptor_base::-4}`
ppkrec=`echo ${receptor_base::-4}'.clean.ppk.pdb'`
echo "DEBUG| " $clean_rec $rec_name $ppkrec
ls ???'_'$rec_name'.pdb' > motif_list
for a in `cat motif_list`; do $MASTER/createPDS --type query --pdb $a; done

n_searches=$(wc -l motif_list | gawk '{print $1}')

# Step 2: Run MASTER
run_master_jid=$(sbatch --array=0-"$n_searches"%50 run_master.sh $master_cutoff | awk '{print $NF}')

# Step2.5: Prepack receptor
prep_input_jid=$(sbatch --dependency=afterany:"${run_master_jid}" --job-name=prep_input --get-user-env --time=90:00:00\
                --mem=1600m prepare_input.sh $clean_rec | awk '{print $NF}')

# Step3: Extract templates
extract_templates_jid=$(sbatch --array=0-"$n_searches"%50 --dependency=afterok:"${prep_input_jid}" run_extract_templates.sh \
                    "$pep_sequence" "$ppkrec" | awk '{print $NF}')

# Step 4: FPD
fpd_jid=$(sbatch --dependency=afterany:"${extract_templates_jid}" --chdir=$(pwd) --job-name=fpd \
          fpd.sh "$clean_rec" "$min_rec_bb" | awk '{print $NF}')


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
#	Notification script: $([ "$notify_script" ] && echo $notify_script)
	------------------------------------------------

	JOBINFO

}

echo "SLURM_JID $finalize_jid"
