/**
 * @file    mpi_findmotifs.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the functions for solving the findmotif problem in
 *          parallel using a master-worker paradigm.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
// You DO NOT need to change anything in this file.
#ifndef MPI_FINDMOTIFS_H
#define MPI_FINDMOTIFS_H

#include "findmotifs.h"


/**
 * @brief   Solves the motif finding problem for a partially solved problem.
 *          The exhaustive enumeration starts from the value given by
 *          `start_value` and enumerates all values which differ in up to
 *          `d- hamming(start_value, input[0]` bits. This function considers
 *          bit inversions only in the bits [startbitpos, l). This function
 *          returns all values that are at most `d` bits different from
 *          all input values.
 *
 * @param n         The number of input sequences.
 * @param l         The length (in bits) of each input sequence.
 * @param d         The number of bits that are allowed to differ between
 *                  the answers and all input sequences.
 * @param input     An array of the `n` input sequences. Represented as
 *                  64 bit integers.
 * @param startbitpos   Restricts the range of bit inversions to bits
 *                      positioned higher (towards MSB) or equal to than this
 *                      one. This is equal to the `master_depth` parameter
 *                      in the master process. The master process
 *                      enumerates the changes to all bits in the range [0,
 *                      startbitpos). Then this function continues in the range
 *                      [startbitpos, l).
 * @param start_value   Instead of starting with `input[0]`, start with this
 *                      value and mutate only `d-hamming(start_value, input[0]`
 *                      further bits.
 *
 * @return          A std::vector containing all answers for this subproblem.
 */
std::vector<bits_t> findmotifs_worker(unsigned int n, unsigned int l,
                                     unsigned int d, const bits_t* input,
                                     unsigned int startbitpos,
                                     bits_t start_value);

/**
 * @brief Solves the motif finding problem till a specified depth
 *        `master_depth`. This function enumerates all values
 *        that are at most `d` bits different from `input[0]` by inverting
 *        bits in the range [0, master_depth) and for each value found,
 *        it passes the remaining problem to a worker process (which mutates
 *        all possible combinations of bits in the range [master_depth, l).
 *
 * @param n         The number of input sequences.
 * @param l         The length (in bits) of each input sequence.
 * @param d         The number of bits that are allowed to differ between
 *                  the answers and all input sequences.
 * @param input     An array of the `n` input sequences. Represented as
 *                  64 bit integers.
 * @param master_depth  The "recursion" depth, until which the master process
 *                      solves the problem before passing the remaining
 *                      subproblem to a worker process. (Note: this function
 *                      does not necessarily have to be solved using
 *                      actual recursion).
 *
 * @return          A std::vector containing all answers.
 */
std::vector<bits_t> findmotifs_master(unsigned int n, unsigned int l,
                                      unsigned int d, const bits_t* input,
                                      unsigned int master_depth);

/**
 * @brief The main worker function. This receives the parameters (n,l,d) and
 *        input parameters from the master process. Then this function waits
 *        for subproblems send by the master process and solves the received
 *        subproblems (using the `findmotifs_worker(...)` function). The
 *        solutions are send back while already solving the next subproblem
 *        (overlapping communication and computation.
 */
void worker_main();


/**
 * @brief   The master process' main function. It solves the complete problem
 *          by communicating with the slave processes and using
 *          the `findmotifs_master(..)` function.
 *
 * @param n         The number of input sequences.
 * @param l         The length (in bits) of each input sequence.
 * @param d         The number of bits that are allowed to differ between
 *                  the answers and all input sequences.
 * @param input     An array of the `n` input sequences. Represented as
 *                  64 bit integers.
 * @param master_depth  The "recursion" depth, until which the master process
 *                      solves the problem before passing the remaining
 *                      subproblem to a worker process. (Note: this function
 *                      does not necessarily have to be solved using
 *                      actual recursion).
 *
 * @return          A std::vector containing all answers.
 */
std::vector<bits_t> master_main(unsigned int n, unsigned int l, unsigned int d,
                                const bits_t* input, unsigned int master_depth);

#endif // MPI_FINDMOTIFS_H
