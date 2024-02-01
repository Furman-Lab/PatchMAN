#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=93:00:00
#SBATCH --mem=1G
#SBATCH --nice=8000
#SBATCH --get-user-env

master_list=($(ls *pds))
echo "$SLURM_ARRAY_TASK_ID"
master_job=${master_list["$SLURM_ARRAY_TASK_ID"]}

# For benchmarking, filter out those pdbs that are in the list of similar PDBs
if [[ $# -eq 2 ]]
then
	targetList=$2
else
	targetList="$PROTOCOL_ROOT/db_list_30nonred"
fi

command="$MASTER/master --query $master_job --targetList $targetList --bbRMSD --rmsdCut $1 --topN 100000 --matchOut $master_job.matches"

echo $command

$command

