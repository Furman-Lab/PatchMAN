#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=13:00:00
#SBATCH --mem=1G
#SBATCH --get-user-env

echo ${PYTHON} ${PROTOCOL_ROOT}/bin/prepack_receptor.py $1
${PYTHON} ${PROTOCOL_ROOT}/bin/prepack_receptor.py $1
