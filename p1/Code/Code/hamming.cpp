// Implement your solutions in this file
#include "hamming.h"

// Lookup Table for a 8-bit number. 
// BitsSetTable256[n] = # of 1 in bit representation of n; for n <= 255.
static const unsigned char BitsSetTable256[256] = 
{
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
        B6(0), B6(1), B6(1), B6(2)
};


unsigned int hamming(uint64_t x, uint64_t y)
{
    // TODO: calculate the hamming distance and return it
    
    // XOR operation x = x^y
    x ^= y;

    return BitsSetTable256[x&0xff] + BitsSetTable256[(x>>8)&0xff] + BitsSetTable256[(x>>16)&0xff] +
    BitsSetTable256[(x>>24)&0xff] + BitsSetTable256[(x>>32)&0xff] + BitsSetTable256[(x>>40)&0xff] +
    BitsSetTable256[(x>>48)&0xff] + BitsSetTable256[(x>>56)&0xff];
}
