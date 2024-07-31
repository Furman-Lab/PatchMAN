#!/bin/bash
#SBATCH --job-name=extract_temp
#SBATCH --time=93:00:00
#SBATCH --mem=1600m
#SBATCH --nice=8000
#SBATCH --get-user-env

matches=($(ls *matches))

match_list=${matches["$SLURM_ARRAY_TASK_ID"]}
motifs=($(cat motif_list))
motif=${motifs["$SLURM_ARRAY_TASK_ID"]}

echo $1 > pepfile

# For old DB uncomment:
echo ${PYTHON}/extract_peps_for_motif.py -m "$match_list" -p pepfile -r "$2" --patch "$motif"
${PYTHON} ${BIN_DIR}/extract_peps_for_motif.py -m "$match_list" -p pepfile -r "$2" --patch "$motif"

