#!/usr/bin/make -f

CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline

all: example_floatbox example_velocityboxfiler

example_floatbox: example_floatbox.c
	$(CC) $(CCFLAGS) -o example_floatbox example_floatbox.c

example_velocityboxfiler: example_velocityboxfiler.c
	$(CC) $(CCFLAGS) -o example_velocityboxfiler example_velocityboxfiler.c

clean:
	$(RM) -f example_floatbox example_velocityboxfiler
