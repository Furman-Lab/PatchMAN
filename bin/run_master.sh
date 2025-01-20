#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=93:00:00
#SBATCH --mem=1G
#SBATCH --nice=8000
#SBATCH --get-user-env

master_list=($(ls *pds))
echo "$SLURM_ARRAY_TASK_ID"
master_job=${master_list["$SLURM_ARRAY_TASK_ID"]}

master --query $master_job --bbRMSD --rmsdCut $1 --targetList $2 --topN 100000 --matchOut $master_job.matches


