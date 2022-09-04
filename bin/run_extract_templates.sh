#!/bin/bash
#SBATCH --job-name=extract_temp
#SBATCH --time=93:00:00
#SBATCH --mem=1600m
#SBATCH --nice=8000

matches=($(ls *matches))

match_list=${matches["$SLURM_ARRAY_TASK_ID"]}
motifs=($(cat motif_list))
motif=${motifs["$SLURM_ARRAY_TASK_ID"]}

echo $1 > pepfile

python ${BIN_DIR}/extract_templates.py -m "$match_list" -p pepfile -r "$2" --patch "$motif"

#echo 'TEST'