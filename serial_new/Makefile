CC ?= cc
CFLAGS = -O3 -Wfatal-errors
DOCS = ../docs
RM ?= rm

all: sweep-tt-multistart

sweep-tt-multistart: sweep-tt-multistart.c
	$(CC) -o sweep-tt-multistart $(CFLAGS) sweep-tt-multistart.c

# starting point is 0 based so subtract 1
run: sweep-tt-multistart
	./sweep-tt-multistart $(DOCS)/velocity-241-241-51.txt $(DOCS)/818-FS.txt $(DOCS)/start-4-241-241-51.txt

clean:
	$(RM) -f *.o sweep-tt-multistart
