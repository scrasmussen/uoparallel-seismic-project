#!/usr/bin/make -f

CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline

all: example_floatbox example_iovelocity

example_floatbox: example_floatbox.c
	$(CC) $(CCFLAGS) -o example_floatbox example_floatbox.c

example_iovelocity: example_iovelocity.c
	$(CC) $(CCFLAGS) -o example_iovelocity example_iovelocity.c

clean:
	$(RM) -f example_floatbox example_iovelocity
