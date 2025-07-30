#!/bin/bash
#SBATCH	--mem-per-cpu=2gb
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nice=8000
#SBATCH	--kill-on-invalid-dep=yes
#SBATCH	--get-user-env
#SBATCH --job-name=fpd

echo "FPD called with arguments: $@"
time python3 "$PROTOCOL_ROOT/bin/fpd.py" "$@" 2>&1 | grep -vi WARNING #| grep -v 'Rosetta thread index was requested'
