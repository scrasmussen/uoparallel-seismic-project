===IMPORTANT=== 
If the current directory is somewhere other than $HOME/uoparallel-seismic-project/mpi you will need
to edit the DIR varialbe in run.sh to point where the mpi directory is.

run.sh : this file allocates compute node from aciss and run it by mpirun.This file does not need to be executed.

The main code is in the file sweep-tt-multistart.c file. In the code, there are two matrices for dividing the the whole data(sizeOfTask is the array name) and there is a communication matrix to define the send and receive tasks and ghost regions. 

There is a spreadsheet called mpi.ods which contains how the communication matrix was generated.

There are two files called acissparallel and acissserial which contains the output of running the program serially and parallely.
