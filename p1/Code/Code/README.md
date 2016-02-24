CSE 6220 Programming Assignment 1
=================================

The PDF file supplied with the assignment explains the problem you have to solve
using a master-worker paradigm.

We provide this framework to get you started. There is a `main.cpp` file, which
is the main entry point to the program you have to write. This is already
fully implemented. This function reads command line parameters, reads the
input file and finally writes the output file. You should take a look at this
file, but you won't have to change anything in it. We also provide a program
`generate_input.cpp` that can generate random input for the problem at hand.

For this programming assignment, you will have to implement the functions
that are in the files `hamming.cpp`, `findmotifs.cpp`, and `mpi_findmotifs.cpp`.
The functions that you have to implement are declared and documented in the
corresponding `.h` files.


Compiling
---------
We provide a makefile for compiling the code. Simply run `make` in the
projects directory to build the two executables: `findmotifs` and
`generate_input`. Run each without commandline parameters for seeing usage
descriptions.

Running on Jinx
---------------
To run your parallel program on jinx, we supply a very basic PBS script
in `pbs_script.sh`. You'll have to change this script to point to the correct
project folder. You'll also have to modify the script to achieve multiple runs
for your experiments and for testing different inputs. DO NOT just take this
script as it is.

Sample Input and Output
-----------------------
You can find sample input and the corresponding solution files in the `input`
folder.  For this project, you will have to also solve different inputs
generated the `generate_input` function. Try to find some settings of
parameters, such that your implementation takes between 5-10 minutes
sequentially and show your speedups for such an instance. Make sure the output
does not become too large.  Generally speaking, the size of the output decreases
with increasing `n`, while it drastically increases with increasing `d`. Play
with all three parameters (n,l,d) to find instances that have a long enough run
time but still have small enough output (approx < 1 MiB).
