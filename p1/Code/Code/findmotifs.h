/**
 * @file    findmotifs.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the sequential findmotifs function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
// You DO NOT need to change anything in this file.
#ifndef FINDMOTIFS_H
#define FINDMOTIFS_H

#include <vector>
#include <stdint.h>

/// The datatype used for representing a (up to) 64 bit sequence.
typedef uint64_t bits_t;

/**
 * @brief   Solves the motif finding problem _sequentially_.
 *
 * @param n         The number of input sequences.
 * @param l         The length (in bits) of each input sequence.
 * @param d         The number of bits that are allowed to differ between
 *                  the answers and all input sequences.
 * @param input     An array of the `n` input sequences. Represented as
 *                  64 bit integers.
 *
 * @return          A std::vector containing all answers.
 */
std::vector<bits_t> findmotifs(unsigned int n, unsigned int l, unsigned int d,
                               const bits_t* input);

/**
 * @brief  Consider bit-representation of an unsigned int A as the modification scheme, 
 *         which shows at which bit of a given sequence S is changed. After the modification, 
 *                                 Result = S XOR A
 *
 *         Given a existed modification position scheme A, this function is going to give the next 
 *         modificaiton scheme following the lexicographic order. 
 *
 * @param scheme    Modification scheme, recorded as an unsigned int.
 * @param ones     The number of 1 in the bit of given 'scheme'.
 * @param seq_len   The length of scheme, only last seq_len bits of the scheme are used
 * @param bound     The upperbound of the number of one in the scheme
 *
 * @return          The number of 1 in the next scheme. New scheme is directly stored in the 'scheme' 
 *                  If 'scheme' is the last one, return negative number.
 */

int nextModiScheme(bits_t &scheme,  
                     unsigned int ones, unsigned seq_len,
                     unsigned int bound);

#endif // FINDMOTIFS_H
