/**
 * @file    io.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implementation of common IO routines.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "io.h"

std::vector<int> scatter_stream_block_decomp_slave(MPI_Comm comm)
{
    // get MPI Communicator properties
    int rank, p;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    // the local vector size (MPI restricts message sizes to `int`)
    int local_size;

    // init result
    std::vector<int> local_elements;

    if (rank == 0)
    {
        std::runtime_error("master called slave function");
    }
    else
    {
        /* I am NOT the root process */

        unsigned long long n;
        // bcast the global size
        MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
        local_size = block_decompose(n, p, rank);

        // resize local data
        local_elements.resize(local_size);
        // actually receive the data
        MPI_Recv (&local_elements[0], local_size, MPI_INT, 0, rank, comm, MPI_STATUS_IGNORE);
    }

    // return the local vectors
    return local_elements;
}

std::vector<int> scatter_file(const char* filepath, MPI_Comm comm)
{
    // init result
    std::vector<int> local_elements;

    // get rank
    int rank;
    MPI_Comm_rank(comm, &rank);

    if (rank == 0)
    {
        /* I am the root process */
        // get file stream iterator
        std::ifstream infile_stream(filepath, std::ifstream::in);
        if (!(infile_stream.good() && infile_stream.is_open()))
                throw std::runtime_error(std::string("Couldn't open file ") + filepath);
        std::istream_iterator<int> input_iterator(infile_stream);

        // get number of elements as first element from stream
        int n = *(input_iterator++);

        // scatter the file accross processors
        local_elements = scatter_stream_block_decomp(input_iterator, n, comm);
    }
    else
    {
        /* I am NOT the root process */
        local_elements = scatter_stream_block_decomp_slave(comm);
    }

    return local_elements;
}

std::vector<int> scatter_vector_block_decomp(std::vector<int>& global_vec, MPI_Comm comm)
{
    // get MPI Communicator properties
    int rank, p;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    // the local vector size (MPI restricts message sizes to `int`)
    int local_size;

    // init result
    std::vector<int> local_elements;

    if (rank == 0)
    {
        /* I am the root process */

        // get size of global array and bcast
        unsigned long long n = global_vec.size();
        MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
        local_size = block_decompose(n, p, rank);

        // scatter-v the actual data
        std::vector<int> counts(p);
        std::vector<int> displs(p);
        displs[0] = 0;
        for (int i = 0; i < p; ++i)
        {
            counts[i] = block_decompose(n, p, i);
            if (i > 0)
                displs[i] = displs[i-1] + counts[i-1];
        }
        local_elements.resize(local_size);
        MPI_Scatterv(&global_vec[0], &counts[0], &displs[0], MPI_INT,
                     &local_elements[0], local_size, MPI_INT, 0, comm);
    }
    else
    {
        /* I am NOT the root process */

        // receive the size of my local array
        unsigned long long n;
        MPI_Bcast(&n, 1, MPI_UNSIGNED_LONG_LONG, 0, comm);
        local_size = block_decompose(n, p, rank);

        // resize result buffer
        local_elements.resize(local_size);
        // actually receive all the data
        MPI_Scatterv(NULL, NULL, NULL, MPI_INT,
                     &local_elements[0], local_size, MPI_INT, 0, comm);
    }

    // return local array
    return local_elements;
}

std::vector<int> gather_vectors(std::vector<int>& local_vec, MPI_Comm comm)
{
    // get MPI parameters
    int rank;
    int p;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &p);

    // get local size
    int local_size = local_vec.size();

    // init result
    std::vector<int> result;

    // master process: receive results
    if (rank == 0)
    {
        // gather local array sizes, sizes are restricted to `int` by MPI anyway
        // therefore use int
        std::vector<int> local_sizes(p);
        MPI_Gather(&local_size, 1, MPI_INT,
                   &local_sizes[0], 1, MPI_INT,
                   0, comm);

        // gather-v to collect all the elements
        int total_size = std::accumulate(local_sizes.begin(), local_sizes.end(), 0);
        result.resize(total_size);

        // get receive displacements
        std::vector<int> displs(p, 0);
        for (int i = 1; i < p; ++i)
            displs[i] = displs[i-1] + local_sizes[i-1];

        // gather v the vector data to the root
        MPI_Gatherv(&local_vec[0], local_size, MPI_INT,
                    &result[0], &local_sizes[0], &displs[0], MPI_INT, 0, comm);
    }
    // else: send results
    else
    {
        // gather local array sizes
        MPI_Gather(&local_size, 1, MPI_INT, NULL, 1, MPI_INT, 0, comm);

        // sent the actual data
        MPI_Gatherv(&local_vec[0], local_size, MPI_INT,
                    NULL, NULL, NULL, MPI_INT, 0, comm);
    }

    return result;
}
