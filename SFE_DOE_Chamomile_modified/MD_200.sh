#!/bin/bash 
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --cpus-per-task=32
#SBATCH --output=matlab_parallel_200_n.out

module load matlab

srun matlab -nodisplay -r 'MD(200, 32, 32)'
