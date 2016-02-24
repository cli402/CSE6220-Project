/**
 * @file    utils.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements common utility/helper functions.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *             You can add new functions to this header.             *
 *********************************************************************/

#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

/*********************************************************************
 * DO NOT CHANGE THE FUNCTION SIGNATURE OF THE FOLLOWING 2 FUNCTIONS *
 *********************************************************************/

inline size_t block_decompose(const size_t n, const int p, const int rank)
{
    return n / p + (((size_t)rank < n % p) ? 1 : 0);
}

inline size_t block_decompose(const size_t n, MPI_Comm comm)
{
    int rank, p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);
    return block_decompose(n, p, rank);
}

/*********************************************************************
 *                  DECLARE YOUR OWN FUNCTIONS HERE                  *
 *********************************************************************/

// ...


#endif // UTILS_H
