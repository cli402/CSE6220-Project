/**
 * @file    io.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements common IO and input generation functions.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

#ifndef IO_H
#define IO_H

#include <mpi.h>

#include <fstream>
#include <iterator>
#include <vector>
#include <string>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <algorithm>
#include <numeric>

#include "utils.h"


template <typename Iterator>
void fill_with_rand(Iterator in, std::size_t n)
{
    std::generate_n(in, n, rand);
}

template <typename InputIterator, typename OutputIterator>
void copy_n(InputIterator& in, std::size_t n, OutputIterator out)
{
    for (std::size_t i = 0u; i < n; ++i)
        *(out++) = *(in++);
}

template <typename Iterator>
std::vector<int>
scatter_stream_block_decomp(Iterator input, unsigned long long n, MPI_Comm comm)
{
    // get MPI Communicator properties
    int rank, p;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    typedef typename std::iterator_traits<Iterator>::value_type val_t;


    // init result
    std::vector<val_t> local_elements;

    if (rank == 0)
    {
        /* I am the root process */
        // the local vector size (MPI restricts message sizes to `int`)
        unsigned long long local_size_l = block_decompose(n, p, rank);
        if (local_size_l > std::numeric_limits<int>::max())
            throw std::runtime_error("Local size must be smaller than int limit");
        int local_size = local_size_l;

        // bcast the global size
        MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);

        // copy the first block into the masters memory
        local_elements.resize(local_size);
        copy_n(input, local_size, local_elements.begin());

        // distribute the data one-by-one
        // (since we don't have a parallel filesystem)
        for (int i = 1; i < p; ++i) {
            // copy into local buffer
            std::vector<val_t> local_buffer(local_size);
            int remote_size = block_decompose(n, p, i);
            copy_n(input, remote_size, local_buffer.begin());
            // send the data to processor i
            MPI_Send (&local_buffer[0], remote_size, MPI_INT, i, i, comm);
        }
    }
    else
    {
        /* I am NOT the root process */
        std::runtime_error("slave called master function");
    }

    // return the local vectors
    return local_elements;
}

std::vector<int> scatter_stream_block_decomp_slave(MPI_Comm comm);

std::vector<int> scatter_file(const char* filepath, MPI_Comm comm);

std::vector<int> scatter_vector_block_decomp(std::vector<int>& global_vec, MPI_Comm comm);

std::vector<int> gather_vectors(std::vector<int>& local_vec, MPI_Comm comm);

#endif // IO_H
