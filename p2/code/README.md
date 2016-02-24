CSE 6220 Programming Assignment 2
=================================

## Code hosting

You are highly encouraged to use a code version management tool such as git.
This will help you to code in a team and keep track of your progress.

However, do not upload your code to public repositories. If you want to use
version management for collaboration, make sure to use private repositories.

I highly recommend using Georgia Tech's Enterprise Github installation at
https://github.gatech.edu/

We will also be hosting the framework code for the programming assignment on
there.  If you find any issues with the code framework and we have to make
changes, we will publish those changes in that GitHub repository additionally to
sending out the updated framework.

## Code structure

All the code is located at the root level of the project.
The `gtest` folder contains the Google Test Unit Testing framework version 1.7.

There are multiple header and .cpp files, your implementation will go
into the following files:

- `jacobi.cpp`: Implement the sequential algorithm for Jacobi's method according
  to the function declarations in `jacobi.h`
- `mpi_jacobi.cpp`: Implement the parallel algorithm according to the
  declarations in `mpi_jacobi.h`.
- `utils.h` and `utils.cpp`: Implement common utility functions in these files
- `seq_tests.cpp`: Unit tests for the sequential code, implement your own
  test cases in here
- `mpi_tests.cpp`: Unit tests for the parallel MPI code. Implement your own
  test cases for all declared functions in here.


Other files containing code that you should not change are:

- `main.cpp`: Implements code for the main executable `jacobi`. This does
  input/output reading and calling of the actual functions.
- `io.h`: implements IO functions and random input generation
- `mpi_gtest.cpp`: MPI wrapper for the GTest framework


Utility scripts (you may play around with these to generate your own custom
input):

- `generate_input.py`: Python script to generate inputs. You can modify this
  code to generate different inputs.
- `check_output.py`: Checks whether the output from `jacobi` is correct, by
  comparing the output with pythons numpy implementation.


## Compiling

In order to compile everything, simply run
```sh
make all
```


## Running Tests

For running all tests do:
```sh
make test
```

You can also run the two tests (sequential/MPI) separately by either:
```sh
./seq_tests
```
or
```sh
mpirun -np 4 ./mpi_tests
```
