# Makefile for HPC 6220 Programming Assignment 1
CC=mpic++
#CCFLAGS=-Wall -g
# activate for compiler optimizations:
CCFLAGS=-Wall -O3 -g
LDFLAGS=
 
all: generate_input findmotifs

generate_input: generate_input.o
	$(CC) $(LDFLAGS) -o $@ $^

findmotifs: main.o findmotifs.o mpi_findmotifs.o hamming.o
	$(CC) $(LDFLAGS) -o $@ $^

%.o: %.cpp %.h
	$(CC) $(CCFLAGS) -c $<

%.o: %.cpp
	$(CC) $(CCFLAGS) -c $<

clean:
	rm -f *.o generate_input findmotifs
