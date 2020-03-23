#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem=20GB
#SBATCH --time=30:00:00
#SBATCH --job-name=MTHMTC
#SBATCH --mail-type=ALL

module purge
module load mathematica/11.0.1
srun math < ServerScript.wl

