#!/usr/bin/make -f

CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline
INCLUDE = -I../include
RM ?= rm

all: vconvert

vconvert: vconvert.c
	$(CC) $(CCFLAGS) -o vconvert vconvert.c $(INCLUDE)

clean:
	$(RM) -f vconvert
