Instructions:

I have deleted the velocity file velocity-241-241-51.txt to make the size smaller.Please copy this file before running at aciss.

There are two scripts.
1. all.sh : please run this script in ACISS, it compiles and call the
run.sh file. Then it will store the output in output.tt file and all
response from program will be stored in "output" file.

2. run.sh : this file allocates compute node from aciss and run it by mpirun.This file does not need to be executed. This file will be called by all.sh file.

The main code is in the file sweep-tt-multistart.c file. In the code, there are two matrices for dividing the the whole data(sizeOfTask is the array name) and there is a communication matrix to define the send and receive tasks and ghost regions. 

There is a spreadsheet called mpi.ods which contains how the communication matrix was generated.

There are two files called acissparallel and acissserial which contains the output of running the program serially and parallely.
