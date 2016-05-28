#!/usr/bin/make -f

CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline
INCLUDE = -I../include

all: example_floatbox example_velocityboxfiler example_mpihelpers

example_floatbox: example_floatbox.c
	$(CC) $(CCFLAGS) -o example_floatbox example_floatbox.c $(INCLUDE)

example_velocityboxfiler: example_velocityboxfiler.c
	$(CC) $(CCFLAGS) -o example_velocityboxfiler example_velocityboxfiler.c $(INCLUDE)

example_mpihelpers: example_mpihelpers.c
	$(CC) $(CCFLAGS) -o example_mpihelpers example_mpihelpers.c

clean:
	$(RM) -f example_floatbox example_velocityboxfiler example_mpihelpers
