CC = mpicc
CFLAGS = -O3
DOCS=../docs

all: sweep-tt-multistart

sweep-tt-multistart: sweep-tt-multistart.c
	$(CC) -o sweep-tt-multistart $(CFLAGS) sweep-tt-multistart.c

job: sweep-tt-multistart
	qsub -q generic run.sh

clean:
	rm -f run.sh.* output sweep-tt-multistart output.tt *.o
