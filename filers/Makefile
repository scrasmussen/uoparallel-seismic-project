C ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline

all: iotest

iotest: iotest.c iovelocity.h
	$(CC) $(CCFLAGS) -o iotest iotest.c

clean:
	$(RM) -f iotest
