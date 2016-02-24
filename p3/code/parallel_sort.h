/**
 * @file    parallel_sort.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Declares the parallel sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#ifndef PARALLEL_SORT_H
#define PARALLEL_SORT_H

#include <mpi.h>

/**
 * @brief   Parallel, distributed sorting over all processors in `comm`. Each
 *          processor has the local input [begin, end).
 *
 * Note that `end` is given one element beyond the input. This corresponds to
 * the API of C++ std::sort! You can get the size of the local input with:
 * int local_size = end - begin;
 *
 * @param begin Pointer to the first element in the input sequence.
 * @param end   Pointer to one element past the input sequence. Don't access this!
 * @param comm  The MPI communicator with the processors participating in the
 *              sorting.
 */
void parallel_sort(int* begin, int* end, MPI_Comm comm);




/**
 * @brief   Parallel, distributed sorting over all processors in `comm`. Each
 *          processor has the local input [begin, begin+arraysize). Besizdes,
 *          the size of the local array may change along the solving. The new
 *          size is returned by the function
 *
 * Note that `end` is given one element beyond the input. This corresponds to
 * the API of C++ std::sort! You can get the size of the local input with:
 * int local_size = end - begin;
 *
 * @param begin Pointer to the first element in the input sequence.
 * @param size  integer showing the present length of the array!
 * @param comm  The MPI communicator with the processors participating in the
 *              sorting.
 *
 * @return The new size of local array
 */

int parallel_sort_changeable(int* &begin, int size, MPI_Comm comm);


/**
 * @brief   Partition the processors for the two separated array
 *
 * @param procLeft   integer, used to get the number of processor for the left-part array
 * @param procRight  integer, used to get the number of processor for the right-part array
 * @param numProc    integer, total number of processor  
 * @param numLeft    integer, the number of smaller-than-or-equal-to-pivot elements
 * @param numRight   integer, the number of larger-than-pivot elements
 */
void Partition_Proc(int &procLeft, int &procRight, int numProc, int numLeft, int numRight, MPI_Comm comm);

/**
 * @brief   At the end of sorting, Calcluate the configuration for reallocate the array
 *          for present processor based on the new_localsize and object local size;
 *          Besides, we assume that the total amount of elements remain the same based on 
 *          the given presentSize and ObjectiveSize in each processor.
 *
 * @param sendcnts     Pointer, Send Counts, array with length numProc, 
 *                     i th entry stores the number of elements received from Processor i; 
 * @param sdispls      Pointer, Send Displacement, array with length numProc,
 *                     i th entry stores the displacement of the sendbuf for Processor i;
 * @param recvcnts     Pointer, Receive Counts, Similar to sendcnts;
 * @param rdispls      Pointer, Receive Displacement, Similar to sdispls;
 * @param curSize      Integer, current size of the local array
 * @param objectSize   Integer, the required final local array size at each processor
 * @param comm         MPI_Comm, the communicator involved in this computation
 */

void Assign_Data(int* sendcnts, int* sdispls, 
            int* recvcnts, int* rdispls, int curSize, int objectSize, MPI_Comm &comm);

/**
 * @brief   Partition the local array according to the given pivot and store the size of
 *          Left part;
 * @param  begin     pointer, the beginning of the local array
 * @param  end       pointer, the end of the local array;
 * @param  pivot     integer, Pivot
 * @param  numLeft   integer, used to get the size of left part after partition
 */
void Partition_Local(int* begin, int* end, int pivot, int &numLeft, int &numEqual);


/*********************************************************************
 *              Declare your own helper functions here               *
 *********************************************************************/

// ...

#endif // PARALLEL_SORT_H
