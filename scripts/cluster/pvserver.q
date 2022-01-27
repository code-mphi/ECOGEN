#!/bin/bash
 
#SBATCH --job-name=paraview
#SBATCH -o parav_now.out
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=120000 
#SBATCH --export=ALL 
#SBATCH -t 24:00:00

mpirun /home/kevinsch/software/ParaView-5.6.0-MPI-Linux-64bit/bin/pvserver --force-offscreen-rendering
