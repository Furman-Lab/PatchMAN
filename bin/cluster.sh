#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G
#SBATCH --nodes=1
#SBATCH --get-user-env

mkdir -p clustering/
python3 "$PROTOCOL_ROOT/bin/cluster.py" -c ${SLURM_CPUS_PER_TASK:-1} > clustering/clog 2>&1
