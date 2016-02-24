/**
 * @file    generate_input.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   A program for generating random input for the findmotif
 *          assignment.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
// You DO NOT need to change anything in this file.
#include <cstdlib>
#include <stdint.h>

#include <iostream>
#include <sstream>
#include <fstream>

typedef uint64_t bits_t;

/**
 * @brief Generates a random unsigned 64 bit integer.
 */
uint64_t rand64()
{
    // concatenate the lower significant bits of consecutive calls to rand()
    // use only the lowest 16 bits, since not all compiles will guarantuee
    // RAND_MAX to be larger (~30 bits)
    uint64_t rnd = 0;
    rnd |= ((uint64_t)(rand() & 0xffff) << 48);
    rnd |= ((uint64_t)(rand() & 0xffff) << 32);
    rnd |= ((uint64_t)(rand() & 0xffff) << 16);
    rnd |= (uint64_t)(rand() & 0xffff);
    return rnd;
}

/**
 * @brief Returns a random instance for the find-motif problem given the
 *        size parameters.
 *
 * @param n     The number of sequences to generate.
 * @param l     The number of bits per sequence.
 * @param d     The maximum hamming distance to a common value.
 *
 * @return      An array of `n` input sequences, represented as 64 bit integers
 *              each.
 */
bits_t* getrandom(unsigned int n, unsigned int l, unsigned int d)
{
    bits_t* input = new bits_t[n];
    bits_t center = rand64() % ((bits_t)0x1 << l);
    for (unsigned int i = 0; i < n; ++i)
    {
        bits_t bits = 0;
        // flip `d` random bits
        for (unsigned int j = 0; j < d; ++j)
        {
            bits |= ((bits_t)0x1 << (rand() % l));
        }
        input[i] = bits ^ center;
    }
    return input;
}

/**
 * @brief Prints the usage of the program to stderr.
 */
void printUsage()
{
    std::cerr << "Usage:  ./generate_input n l d <file>" << std::endl;
    std::cerr << "          n       The number of bit sequences generated (>= 2)." << std::endl;
    std::cerr << "          l       The number of bits in each sequence (<= 64)." << std::endl;
    std::cerr << "          d       The max number of changed bits for the output sequences (<= l)." << std::endl;
    std::cerr << "          <file>  The output filename." << std::endl;
    std::cerr << "          seed    (optional) a seed for the random engine. (default=0)" << std::endl;
}

/**
 * @brief The main function: reads the commandline arguments, opens the file
 *        and writes the results.
 */
int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        printUsage();
        exit(EXIT_FAILURE);
    }

    // parse arguments
    unsigned int n, l, d;
    std::istringstream iss;
    iss.str(argv[1]);
    if (!(iss >> n)) {
        std::cerr << "ERROR: Argument n=" << argv[1] << " is not a number." << std::endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    iss.clear();
    iss.str(argv[2]);
    if (!(iss >> l)) {
        std::cerr << "ERROR: Argument l=" << argv[2] << " is not a number." << std::endl;
        printUsage();
        exit(EXIT_FAILURE);
    }
    iss.clear();
    iss.str(argv[3]);
    if (!(iss >> d)) {
        std::cerr << "ERROR: Argument d=" << argv[3] << " is not a number." << std::endl;
        printUsage();
        exit(EXIT_FAILURE);
    }

    // check for valid values
    if (l > 64 || d > l || n < 2)
    {
        std::cerr << "ERROR: Invalid input values, check constraints." << std::endl;
        printUsage();
        exit(EXIT_FAILURE);
    }

    // open output file
    std::ostream* ostr = &std::cout;
    std::ofstream outfile;
    if (argc >= 5)
    {
        outfile.open(argv[4]);
        if (!(outfile.good() && outfile.is_open()))
        {
            std::cerr << "ERROR: Couldn't open file '" << argv[4] << "' for writing." << std::endl;
            printUsage();
            exit(EXIT_FAILURE);
        }
        ostr = &outfile;
    }
    std::ostream& outstream = *ostr;

    // check if the `seed` argument is set
    int seed = 0;
    if (argc == 6)
    {
        iss.clear();
        iss.str(argv[5]);
        if (!(iss >> seed))
        {
            std::cerr << "ERROR: Argument seed=" << argv[5] << " is not a number." << std::endl;
            printUsage();
            exit(EXIT_FAILURE);
        }
    }

    // seed the random number engine (default = 0)
    srand(seed);

    // generate problemset
    bits_t* input = getrandom(n, l, d);

    // write output file
    outstream << n << " " << l << " " << d << std::endl;
    for (unsigned int i = 0; i < n; ++i)
    {
        outstream << input[i] << std::endl;
    }

    // cleanup
    delete[] input;

    return 0;
}
