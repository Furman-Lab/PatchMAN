#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --mem=2000m
#SBATCH --get-user-env

mkdir -p clustering/
python3 "$PROTOCOL_ROOT/bin/cluster.py" > clustering/clog
