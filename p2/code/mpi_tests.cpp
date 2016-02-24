/**
 * @file    mpi_tests.cpp
 * @ingroup group
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   GTest Unit Tests for the parallel MPI code.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
/*
 * Add your own test cases here. We will test your final submission using
 * a more extensive tests suite. Make sure your code works for many different
 * input cases.
 *
 * Note:
 * The google test framework is configured, such that
 * only errors from the processor with rank = 0 are shown.
 */

#include <mpi.h>
#include <gtest/gtest.h>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include "jacobi.h"
#include "mpi_jacobi.h"
#include "utils.h"
#include "io.h"

/**
 * @brief Creates and returns the square 2d grid communicator for MPI_COMM_WORLD
 */
void get_grid_comm(MPI_Comm* grid_comm)
{
    // get comm size and rank
    int rank, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int q = (int)sqrt(p);
    ASSERT_EQ(q*q, p) << "Number of processors must be a perfect square.";

    // split into grid communicator
    int dims[2] = {q, q};
    int periods[2] = {0, 0};
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, grid_comm);
}

// test parallel MPI matrix vector multiplication
TEST(MpiTest, MatrixVectorMult1)
{
    // simple 4 by 4 input matrix
    double A[4*4] = {10., -1., 2., 0.,
        -1., 11., -1., 3.,
        2., -1., 10., -1.,
        0.0, 3., -1., 8.};
    double x[4] =  {6., 25., -11., 15.};
    double y[4];
    double expected_y[4] = {13.,  325., -138.,  206.};
    int n = 4;

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    // testing sequential matrix multiplication
    mpi_matrix_vector_mult(n, A, x, y, grid_comm);

    // checking if all values are correct (up to some error value)
    for (int i = 0; i < n; ++i)
    {
        EXPECT_NEAR(expected_y[i], y[i], 1e-10) << " element y[" << i << "] is wrong";
    }
}


// test parallel MPI matrix vector multiplication
TEST(MpiTest, Jacobi1)
{
    // simple 4 by 4 input matrix
    double A[4*4] = {10., -1., 2., 0.,
                           -1., 11., -1., 3.,
                           2., -1., 10., -1.,
                           0.0, 3., -1., 8.};
    double b[4] =  {6., 25., -11., 15.};
    double x[4];
    double expected_x[4] = {1.0,  2.0, -1.0, 1.0};
    int n = 4;

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);

    // testing sequential matrix multiplication
    mpi_jacobi(n, A, b, x, grid_comm);

    // checking if all values are correct (up to some error value)
    for (int i = 0; i < n; ++i)
    {
        EXPECT_NEAR(expected_x[i], x[i], 1e-5) << " element y[" << i << "] is wrong";
    }
}


/**
 * Test the parallel code and compare the results with the sequential code.
 */
TEST(MpiTest, JacobiCrossTest1)
{
    // test random matrixes, test parallel code with sequential solutions
    std::vector<double> A;
    std::vector<double> b;
    std::vector<double> mpi_x;

    // get grid communicator
    MPI_Comm grid_comm;
    get_grid_comm(&grid_comm);
    int rank;
    MPI_Comm_rank(grid_comm, &rank);

    int n = 36;
    // initialize data only on rank 0
    if (rank == 0)
    {
        A = diag_dom_rand(n);
        b = randn(n, 100.0, 50.0);
    }

    // getting sequential results
    std::vector<double> x;
    if (rank == 0)
    {
        x.resize(n);
        jacobi(n, &A[0], &b[0], &x[0]);
    }

    // parallel jacobi
    if (rank == 0)
      mpi_x.resize(n);
    mpi_jacobi(n, &A[0], &b[0], &mpi_x[0], grid_comm);

    if (rank == 0)
    {
        // checking if all values are correct (up to some error value)
        for (int i = 0; i < n; ++i)
        {
            EXPECT_NEAR(x[i], mpi_x[i], 1e-8) << " MPI solution x[" << i << "] differs from sequential result";
        }
    }
}
