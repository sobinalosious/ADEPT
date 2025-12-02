#!/bin/bash

#$ -M yliu57@nd.edu
#$ -m abe
#$ -N P312036_Modulus
#$ -pe smp 32
#$ -q  long

module load lammps

fsync $SGE_STDOUT_PATH &

#mpiexec -n $NSLOTS lmp_mpi < in.elastic
mpirun -n $NSLOTS lmp < in.elastic
