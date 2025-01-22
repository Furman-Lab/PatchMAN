#!/bin/bash
#SBATCH --job-name=prepack_receptor
#SBATCH --time=90:00:00
#SBATCH --mem=1600m
#SBATCH --get-user-env

python3 ${PROTOCOL_ROOT}/bin/prepack_receptor.py $1
