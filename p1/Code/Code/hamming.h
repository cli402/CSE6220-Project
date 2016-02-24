/**
 * @file    hamming.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the function calculating the Hamming distance.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
// You DO NOT need to change anything in this file.
#include <stdint.h>
#include <limits.h>

#ifndef HAMMING_H
#define HAMMING_H

/**
 * @brief Returns the hamming distance between x and y.
 *
 * @return  The number of bits that differ between x and y.
 */
unsigned int hamming(uint64_t x, uint64_t y);


#endif // HAMMING_H
