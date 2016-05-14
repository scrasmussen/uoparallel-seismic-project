CC ?= cc
CFLAGS = -O3 -Wfatal-errors
DOCS = ../docs
DATA = ../data
RM ?= rm

all: sweep-tt-multistart

sweep-tt-multistart: sweep-tt-multistart.c
	$(CC) -o sweep-tt-multistart $(CFLAGS) sweep-tt-multistart.c

# starting point is 0 based so subtract 1
run: sweep-tt-multistart
	./sweep-tt-multistart $(DATA)/velocity-241-241-51.vbox $(DOCS)/818-FS.txt $(DOCS)/start-4-241-241-51.txt

clean:
	$(RM) -f *.o sweep-tt-multistart
