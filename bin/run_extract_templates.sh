#!/bin/bash
#SBATCH --job-name=extract_temp
#SBATCH --time=3:00:00
#SBATCH --mem=1600m
#SBATCH --nice=8000
#SBATCH --get-user-env

#$1: peptide sequence
#$2: receptor pdb
#$3: other arguments passed by main script

#source $PROTOCOL_ROOT/.env ",/$(realpath "$work_dir" | cut -d'/' -f2)"

matches=($(ls *matches))

match_list=${matches["$SLURM_ARRAY_TASK_ID"]}
motifs=($(cat motif_list))
motif=${motifs["$SLURM_ARRAY_TASK_ID"]}

echo $1 > pepfile

echo ${PYTHON} ${BIN_DIR}/extract_peps_for_motif.py -m "$match_list" -p pepfile -r "$2" --patch "$motif" $3
${PYTHON} ${BIN_DIR}/extract_peps_for_motif.py -m "$match_list" -p pepfile -r "$2" --patch "$motif" $3
