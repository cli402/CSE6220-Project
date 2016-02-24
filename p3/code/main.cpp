/**
 * @file    main.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the main function which reads input data and dispatches
 *          the sorting.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

/*********************************************************************
 *                  !!  DO NOT CHANGE THIS FILE  !!                  *
 *********************************************************************/

// C std lib
#include <stdint.h>
#include <cstdlib>

// for in/output
#include <iostream>
#include <fstream>
#include <iterator>
#include <vector>
#include <algorithm>

// for MPI
#include <mpi.h>

// own includes
#include "io.h"
#include "utils.h"
#include "parallel_sort.h"

#define DEBUG(msg) std::cerr << msg << std::endl;


/**
 * @brief Returns whether the given input range is in sorted order.
 *
 * @param begin An iterator to the begin of the sequence.
 * @param end   An iterator to the end of the sequence.
 *
 * @return  Whether the input sequence is sorted.
 */
template <typename Iterator>
bool is_sorted(Iterator begin, Iterator end)
{
    if (begin == end)
        return true;
    typename std::iterator_traits<Iterator>::value_type last = *(begin++);
    while (begin != end)
    {
        if (last > *begin)
            return false;
        last = *(begin++);
    }
    return true;
}


/**
 * @brief   Prints a solution.
 *
 * @param solution  The solution to be printed.
 * @param os    The output stream to write to (default: stdout)
 */
void print_solution(std::vector<int>& solution, std::ostream& os = std::cout)
{
    // print out the sorted data
    // this should only run on the master because we do not have parallel IO
    // capabilities on the cluster.
    for (unsigned int i = 0; i < solution.size(); ++i)
    {
        os << solution[i] << std::endl;
    }
}

/**
 * Prints the usage of the program.
 */
void print_usage()
{
    std::cerr << "Usage: ./sort [options] [input_file]" << std::endl;
    std::cerr << "      Optional arguments:" << std::endl;
    std::cerr << "          -o <file>    Output all solutions to the given file." << std::endl;
    std::cerr << "          -o -         Output all solutions to stdout." << std::endl;
    std::cerr << "          -t           Runs global tests on the sorting algorithm. NO input file!" << std::endl;
    std::cerr << "          -r           Run random number tests, random numbers are generated only locally, no bottleneck at startup." << std::endl;
    std::cerr << "          -n <n>       Sets the number of generated input integers (mandatory with option -r and -t)" << std::endl;
    std::cerr << "      Example:" << std::endl;
    std::cerr << "          ./sort -o - input.txt" << std::endl;
    std::cerr << "                  Sorts the numbers given by input.txt and outputs the sorted result to the terminal" << std::endl;
}


int main(int argc, char *argv[])
{
    // Init MPI
    MPI_Init(&argc, &argv);

    /***************************
     *  Parse input arguments  *
     ***************************/

    // forget about first argument (which is the executable's name)
    argc--;
    argv++;

    bool do_output = false;
    char* outfile_path = (char*)"-";

    bool do_read_from_file = true;
    bool do_generate_global = false;
    bool do_generate_local = false;
    int n = -1;
    char* infile_path;

    // parse optional parameters
    while(argc > 0 && argv[0][0] == '-')
    {
        char option = argv[0][1];
        switch (option)
        {
            case 'o':
                // output all solutions to std out
                do_output = true;
                // next parameter is file name
                outfile_path = argv[1];
                argv++;
                argc--;
                break;
            case 't':
                // generate global input for testing
                do_generate_global = true;
                do_read_from_file = false;
                break;
            case 'r':
                if (do_generate_global){
                    print_usage();
                    exit(EXIT_FAILURE);
                }
                do_generate_local = true;
                do_read_from_file = false;
                break;
            case 'n':
                // the next argument must be the number
                argv++;
                argc--;
                n = atoi(argv[0]);
                break;
            default:
                print_usage();
                exit(EXIT_FAILURE);
        }
        // iterate to next argument
        argv++;
        argc--;
    }

    // check that the mandatory parameters are present
    if (do_read_from_file && argc < 1)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }
    if ((do_generate_local || do_generate_global) && n < 1)
    {
        print_usage();
        exit(EXIT_FAILURE);
    }


    // parse mandatory parameters
    if (do_read_from_file)
    {
        infile_path = argv[0];
    }

    /********************
     *  Setting up MPI  *
     ********************/

    // set up MPI
    int rank, p;
    // get total size of processors and current rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);


    // seed the random generator
    // NOTE: This is for random input generation.
    //       During the parallel sort, a fixed seed is used across all processors
    srand(rank);

    /********************************
     *  Distribute input from file  *
     ********************************/
    std::vector<int> local_elements;
    std::vector<int> input_vector;
    if (do_read_from_file)
    {
        // scatter elements from file
         local_elements = scatter_file(infile_path, MPI_COMM_WORLD);
    }
    else if (do_generate_global)
    {
        if (rank == 0)
        {
            // generate huge input array
            input_vector.resize(n);
            fill_with_rand(input_vector.begin(), n);
        }

        // scatter from vector
        local_elements = scatter_vector_block_decomp(input_vector, MPI_COMM_WORLD);
    }
    else if (do_generate_local)
    {
        // fill each local element with randomnes
        int local_size = block_decompose(n, p, rank);
        local_elements.resize(local_size);
        fill_with_rand(local_elements.begin(), local_size);
    }
    else
    {
        print_usage();
        exit(EXIT_FAILURE);
    }

    /**********************************
     *  Parallel, distributed sorting  *
     **********************************/
    // Timing only on master node (and barrierized)
    MPI_Barrier (MPI_COMM_WORLD);
    // start timer
    struct timespec t_start, t_end;
    clock_gettime(CLOCK_MONOTONIC,  &t_start);

    /*****************
     *  Run Sorting  *
     *****************/
    parallel_sort(&local_elements[0], &local_elements[0] + local_elements.size(), MPI_COMM_WORLD);


    MPI_Barrier (MPI_COMM_WORLD);
    // get elapsed time in seconds
    clock_gettime(CLOCK_MONOTONIC,  &t_end);
    double time_secs = (t_end.tv_sec - t_start.tv_sec)
        + (double) (t_end.tv_nsec - t_start.tv_nsec) * 1e-9;
    // output time
    if (rank == 0) {
        DEBUG("sorting took: " << time_secs << " ms");
    }

    /*********************
    *  Collect results  *
    *********************/

    // collect results back to the master node
    std::vector<int> sorted_elements;
    // only gather if output is needed
    if (do_output || do_generate_global)
    {
        sorted_elements = gather_vectors(local_elements, MPI_COMM_WORLD);
    }

    /*****************************
     *  Check or output results  *
     *****************************/

    if (rank == 0)
    {
        if (do_output)
        {
            if (std::string(outfile_path) == "-")
            {
                print_solution(sorted_elements);
            }
            else
            {
                // open file and output to fstream
                std::ofstream outstr(outfile_path);
                print_solution(sorted_elements, outstr);
                outstr.close();
            }
        }

        // if global test, run tests to check correctness
        if (do_generate_global)
        {
            // check that the output is actually sorted
            if (!is_sorted(sorted_elements.begin(), sorted_elements.end()))
            {
                DEBUG("ERROR: OUTPUT IS NOT SORTED");
                exit(EXIT_FAILURE);
            }

            DEBUG("Sorting all sequentially for testing");
            std::sort(input_vector.begin(), input_vector.end());

            DEBUG("Comparing arrays");
            if (! std::equal(input_vector.begin(), input_vector.end(), sorted_elements.begin()))
            {
                DEBUG("ERROR: sorted sequences are not equal");
            }
            else
            {
                DEBUG("This was a triumph! I'm making a note here: HUGE success!");
            }
        }
    }

    // finish up
    MPI_Finalize();
    return 0;
}


