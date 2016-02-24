/**
 * @file    main.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the main routine for the CSE6220 programming assignment 1.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
// You DO NOT need to change anything in this file.
#include <mpi.h>
#include <time.h> // for clock_gettime()

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>

#include "hamming.h"
#include "findmotifs.h"
#include "mpi_findmotifs.h"


void printUsage()
{
    std::cerr << "Usage:  mpirun -np <#processors> ./findmotifs <master-depth> <inputfile> <outputfile>" << std::endl;
    std::cerr << "        <master-depth>            The depth until which the master processor solves the" << std::endl;
    std::cerr << "                                  problem before handing off the the worker threads." << std::endl;
    std::cerr << "        <inputfile>  (optional)   The input file in the format described in the assignment (default: stdin)." << std::endl;
    std::cerr << "        <outputfile> (optional)   The output file, will contain all solutions as integers (default: stdout)." << std::endl;
}


int main(int argc, char *argv[])
{
    // set up MPI
    MPI_Init(&argc, &argv);

    // get communicator size and my rank
    MPI_Comm comm = MPI_COMM_WORLD;
    int p, rank;
    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &rank);

    /* code */
    //testrandom(10, 32, 9);
    bits_t* input;// = getrandom(10, 32, 9);
    if (rank == 0)
    {
        // read input file
        if (argc < 2)
        {
            printUsage();
            exit(EXIT_FAILURE);
        }

        // get the master depth, (i.e. the number of bits till which the
        // master process solves the problem)
        std::string depthstr = argv[1];
        std::istringstream ssargv(depthstr);
        unsigned int master_depth;
        ssargv >> master_depth;

        // open input file (if given, else default to stdin)
        std::istream* instream_ptr = &std::cin;
        std::ifstream infile;
        if (argc >= 3)
        {
            //std::string infilename(argv[2]);
            infile.open(argv[2]);
            if (!(infile.good() && infile.is_open()))
            {
                printUsage();
                exit(EXIT_FAILURE);
            }
            instream_ptr = &infile;
        }
        std::istream& instream = *instream_ptr;

        // open output file (if given, else use stdout)
        std::ostream* outstream_ptr = &std::cout;
        std::ofstream outfile;
        if (argc == 4)
        {
            //std::string outfilename(argv[3]);
            outfile.open(argv[3]);
            if (!(outfile.good() && outfile.is_open()))
            {
                printUsage();
                exit(EXIT_FAILURE);
            }
            outstream_ptr = &outfile;
        }
        std::ostream& outstream = *outstream_ptr;


        // read input from file
        unsigned int n, l, d;
        instream >> n >> l >> d;
        std::vector<bits_t> inputdata(n);
        for (unsigned int i = 0; i < n; ++i)
        {
            instream >> inputdata[i];
        }
        // get data as pointer (for more C-like interface)
        input = &inputdata[0];

        // prepare results
        std::vector<bits_t> results;

        // start timer
        //   we omit the file loading and argument parsing from the runtime
        //   timings, we measure the time needed by the master process
        struct timespec t_start, t_end;
        clock_gettime(CLOCK_MONOTONIC,  &t_start);
        if (p == 1)
        {
            std::cerr << "[WARNING]: Running the sequential solver. Start with mpirun to execute the parallel version." << std::endl;
            // call the sequential solver
            results = findmotifs(n, l, d, input);
        }
        else
        {
            // call the parallel solver function
            results = master_main(n, l, d, input, master_depth);
        }
        // end timer
        clock_gettime(CLOCK_MONOTONIC,  &t_end);
        // time in seconds
        double time_secs = (t_end.tv_sec - t_start.tv_sec)
                         + (double) (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;//???

        // output time
        std::cerr << time_secs << std::endl;

        // write output file in ascending order
        std::sort(results.begin(), results.end());
        for (std::size_t i = 0; i < results.size(); ++i)
        {
            outstream << results[i] << std::endl;
        }
    }
    else
    {
        worker_main();
    }

    // finalize MPI
    MPI_Finalize();
    return 0;
}
