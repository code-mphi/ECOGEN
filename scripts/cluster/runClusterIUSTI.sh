#!/bin/bash
#
#$ -N chocGoutte2
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
#$ -q Douze
#$ -pe ompi_douze 24
#$ -v OMP_NUM_THREADS=1
#
# Name de l'executable
NAME=./ECOGEN
MYARGS=""

. /local/opt/Modules/3.2.9/init/sh
module load openmpi-intel/1.10.1 intel

MYCOMMAND="mpirun -np $NSLOTS ${NAME} ${MYARGS}"

$MYCOMMAND

