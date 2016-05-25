#!/usr/bin/make -f

CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline

all: example_floatbox example_iovelocity example_mpihelpers

example_floatbox: example_floatbox.c
	$(CC) $(CCFLAGS) -o example_floatbox example_floatbox.c

example_iovelocity: example_iovelocity.c
	$(CC) $(CCFLAGS) -o example_iovelocity example_iovelocity.c

example_mpihelpers: example_mpihelpers.c
	$(CC) $(CCFLAGS) -o example_mpihelpers example_mpihelpers.c

clean:
	$(RM) -f example_floatbox example_iovelocity example_mpihelpers
