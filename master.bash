#!/bin/bash
#
#SBATCH -J testi
#SBATCH -t 00:30:00
#SBATCH --mem=30000
#
module load bioconda python-env/3.5.3 && export CONDA_ENVS_PATH=$WRKDIR/DONOTREMOVE/taito-conda-envs && source activate plot2

python3 main.py
