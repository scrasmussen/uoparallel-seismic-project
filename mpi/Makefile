CC = mpicc
CFLAGS = -O3

all: sweep-tt-multistart

sweep-tt-multistart: sweep-tt-multistart.c
	$(CC) -o sweep-tt-multistart $(CFLAGS) sweep-tt-multistart.c

# starting point is 0 based so subtract 1
check: sweep-tt-multistart
	./sweep-tt-multistart velocity-241-241-51.txt 818-FS.txt start-4.txt

clean:
	rm -f *.o sweep-tt-multistart
