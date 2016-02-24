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

#include "io.h"
#include "parallel_sort.h"

/*********************************************************************
 *                   Add your own test cases here.                   *
 *********************************************************************/
// Other test cases can include:
// - all elements are equal
// - elements are randomly picked
// - elements are sorted inversely
// - number of elements is not divisible by the number of processors
// - number of elements is smaller than the number of processors


// test parallel MPI matrix vector multiplication
TEST(MpiTest, Sort10)
{
    int x_in[10] = {4, 7, 5, 1, 0, 2, 9, 3, 8, 6};
    int y_ex[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    std::vector<int> x(x_in, x_in+10);
    std::vector<int> local_x = scatter_vector_block_decomp(x, MPI_COMM_WORLD);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    parallel_sort(&local_x[0], &local_x[0]+local_x.size(), MPI_COMM_WORLD);

    std::vector<int> y = gather_vectors(local_x, MPI_COMM_WORLD);

    if (rank == 0)
        for (int i = 0; i < 10; ++i) {
            EXPECT_EQ(y_ex[i], y[i]);
        }
}

