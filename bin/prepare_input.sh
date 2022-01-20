#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=13:00:00
#SBATCH --mem=1G

python ${BIN_DIR}/prepack_receptor.py $1