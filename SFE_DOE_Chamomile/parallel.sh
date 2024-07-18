#!/bin/bash 
#SBATCH --time=14:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --cpus-per-task=12
#SBATCH --output=matlab_parallel.out

module load matlab

srun matlab -nodisplay -r 'SFE_map'
