#!/bin/bash
#SBATCH -N 2
module load mpi
mpirun -np 4  mpi/shortest-path/sweep-tt-multistart > mpi/shortest-path/output
