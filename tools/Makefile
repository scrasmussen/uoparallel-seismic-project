#!/usr/bin/make -f

CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline
RM ?= rm

all: vconvert

vconvert: vconvert.c
	$(CC) $(CCFLAGS) -o vconvert vconvert.c

clean:
	$(RM) -f vconvert
