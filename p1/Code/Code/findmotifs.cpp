// Implement your solutions in this file
#include "findmotifs.h"
#include "hamming.h"

// implements the sequential findmotifs function
std::vector<bits_t> findmotifs(unsigned int n, unsigned int l,
                               unsigned int d, const bits_t* input)
{
    // If you are not familiar with C++ (using std::vector):
    // For the output (return value) `result`:
    //                  The function asks you to return all values which are
    //                  of a hamming distance `d` from all input values. You
    //                  should return all these values in the return value
    //                  `result`, which is a std::vector.
    //                  For each valid value that you find (i.e., each output
    //                  value) you add it to the output by doing:
    //                      result.push_back(value);
    //                  Note: No other functionality of std::vector is needed.
    // You can get the size of a vector (number of elements) using:
    //                      result.size()

    // create an empty vector
    std::vector<bits_t> result;
    bits_t scheme = 0;
    int diff = 0;

    while( diff >= 0)
    {
        bool is_sol = true;
        bits_t potential_sol = input[0] ^ scheme;

        for(int i = 1; i < n; i++)
        {
            if(hamming(potential_sol, input[i]) > d)
            {
                is_sol = false;
                break;
            }
        }
        if (is_sol)
            result.push_back(potential_sol);
        diff = nextModiScheme(scheme, diff, l, d);
    }

    return result;
}


int nextModiScheme(bits_t &scheme,  
            unsigned int ones, unsigned seq_len,
            unsigned int bound)
{
    //Find the first zero in the bit of scheme
    int pos_zero = 0;
    while (((scheme >> pos_zero) & 1)!= 0 && pos_zero < seq_len)
      pos_zero++;

    if(pos_zero >= seq_len)
    {
        //Given scheme is the biggest one in the lexicographic order
        return -1;
    }
    else if (pos_zero > 0)
    {
        scheme ^= ((1 << (pos_zero+1)) - 1);
        //scheme = ((scheme >> pos_zero) + 1) << pos_zero;
        return ones - pos_zero + 1;
    }
    else
    {
        //Position of the first zero is at first bit
        //Then find the first appeared one. 
        int pos_one = 0;
        while (((scheme >> pos_one)&1)!= 1 && pos_one < seq_len)
          pos_one++;

        if (pos_one >= seq_len)
        {
            //Scheme == 00000000
            scheme = 1;
            return 1;
        }
        else if (ones < bound)
        {
            //Add a '1' at first bit.
            scheme += 1;
            return ones+1;
        }
        else
        {
            bits_t tmp = (scheme >> pos_one);
            int num = nextModiScheme(tmp, ones, seq_len - pos_one, bound);
            scheme = (tmp << pos_one); 
            return num;
        }
    }
}






