#!/bin/bash
#SBATCH --job-name=prepack_receptor
#SBATCH --time=13:00:00
#SBATCH --mem=1G
#SBATCH --get-user-env

${PYTHON} ${PROTOCOL_ROOT}/bin/prepack_receptor.py $1
