CC ?= gcc
CCFLAGS = -std=c11 -O3 -Wall -Wfatal-errors -Wpedantic -Winline

all: boxtest

boxtest: boxtest.c floatbox.h
	$(CC) $(CCFLAGS) -o boxtest boxtest.c

clean:
	$(RM) -f boxtest
