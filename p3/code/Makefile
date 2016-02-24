# Makefile for HPC 6220 Programming Assignment 3
CXX=mpic++
#CCFLAGS=-Wall -g 
# activate for compiler optimizations:
CCFLAGS=-Wall -O3 -g
LDFLAGS=

# set up google test
GTEST_DIR = ./gtest
CCFLAGS += -I. #-I$(GTEST_DIR)


all: sort tests

test: tests
	for p in 3 4 5; do \
	echo "### TESTING WITH $$p PROCESSES ###"; mpirun -np $$p ./mpi_tests ;\
	done

tests: mpi_tests

sort: main.o io.o parallel_sort.o utils.o
	$(CXX) $(LDFLAGS) -o $@ $^

mpi_tests: mpi_tests.o mpi_gtest.o gtest-all.o parallel_sort.o utils.o io.o
	$(CXX) $(LDFLAGS) -o $@ $^

gtest-all.o : $(GTEST_DIR)/gtest-all.cc $(GTEST_DIR)/gtest.h
	$(CXX) $(CCFLAGS) -c $(GTEST_DIR)/gtest-all.cc

%.o: %.cpp %.h
	$(CXX) $(CCFLAGS) -c $<

%.o: %.cpp
	$(CXX) $(CCFLAGS) -c $<

clean:
	rm -f *.o sort seq_tests mpi_tests
