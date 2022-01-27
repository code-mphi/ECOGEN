#!/bin/bash

#SBATCH --job-name=test
#SBATCH -o run.out
#SBATCH -e run.err
#SBATCH --partition=normal
#SBATCH --export=ALL 
#SBATCH -t 24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24

mpirun /home/kevinsch/ECOGEN/ECOGEN_CFPgroupVersion/ECOGEN/ECOGEN

