CC?=gcc
CFLAGS=-O3 -g -fPIC
CFLAGS_KNL=-O3 -g -fPIC -axMIC-AVX512
LDFLAGS=-lm
OPENMP_FLAG=-fopenmp
DOCS=../docs
DELTA=1

all: wave-multistart

wave-multistart: wave-multistart.o
	$(CC) wave-multistart.o -o wave-multistart ${LDFLAGS} ${OPENMP_FLAG}

wave-multistart.o: wave-multistart.c Makefile
	$(CC) wave-multistart.c -o wave-multistart.o -c ${CFLAGS} ${OPENMP_FLAG}

run-wave: wave-multistart
	export OMP_NUM_THREADS=4
	./wave-multistart $(DOCS)/velocity-241-241-51-nonConst.txt $(DOCS)/818-FS.txt $(DOCS)/start-1-241-241-51.txt $(DELTA)

clean:
	rm -f *.o wave-multistart 
