#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=93:00:00
#SBATCH --mem=1G
#SBATCH --nice=8000
#SBATCH --get-user-env

master_list=($(ls *pds))
echo "$SLURM_ARRAY_TASK_ID"
master_job=${master_list["$SLURM_ARRAY_TASK_ID"]}

echo "$MASTER"/master --query "$master_job" --targetList "$PROTOCOL_ROOT"/db_list_30nonred --bbRMSD --rmsdCut 1.5  --topN 100000 --matchOut "$master_job".matches

# For old DB
$MASTER/master --query "$master_job" --targetList "$PROTOCOL_ROOT"/db_list_30nonred --bbRMSD --rmsdCut 1.5  --topN 100000 --matchOut "$master_job".matches
