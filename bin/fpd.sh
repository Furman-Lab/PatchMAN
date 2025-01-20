#!/bin/bash
#SBATCH	--mem-per-cpu=2gb
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --requeue
#SBATCH --killable
#SBATCH --nice=8000
#SBATCH	--time=120:00:00
#SBATCH	--kill-on-invalid-dep=yes
#SBATCH	--get-user-env
#SBATCH --job-name=fpd
# The script didnt work with ntasks=1 and -c 11, so I set ntasks=11 and cpus-per-task=1

echo "FPD called with arguments: $@"
time python3 "$PROTOCOL_ROOT/bin/fpd.py" "$@" 2>&1 | grep -vi WARNING #| grep -v 'Rosetta thread index was requested'
