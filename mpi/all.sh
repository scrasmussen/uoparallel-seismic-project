module load mpi
rm run.sh.* output sweep-tt-multistart output.tt
touch output.tt 
mpicc -o sweep-tt-multistart sweep-tt-multistart.c 
qsub -q generic run.sh
