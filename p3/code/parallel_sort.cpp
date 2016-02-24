/**
 * @file    parallel_sort.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements the parallel, distributed sorting function.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "parallel_sort.h"
#include "utils.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <ctime>


// implementation of your parallel sorting
/*
 * Three conditions: 
 * 1: p > n
 * construct a n-processor communicator out of comm: new_comm
 * run parallel_sort(int *begin, int* end, new_comm)
 *
 * 2: 1 < p <= n;
 *    a: randomly choose the pivot between 0 and end-begin-1;
 *    b: do the partition locally
 *    c: collect the number of each subarray and get the total number m1 and m2; (for each process, Allgather)
 *    d: Assign new subproblems, calculate the destination of sending and receiving.
 *    e: parallel_sort(begin, end, newnewcomm)
 *
 * 3: p = 1;
 * Use available C++ function to sort. No more work to do. 
 */

void parallel_sort(int* begin, int* end, MPI_Comm comm)
{
    //Basic Information
    int numProc;
    int localSize = end - begin;
    MPI_Comm_size(comm, &numProc);

    //Copy [begin, end) to 'cache'
    int* cache = new int[localSize];
    for(int i = 0; i < localSize; i++)
      cache[i] = *(begin+i);

    std::srand(time(NULL));
    /*
    int myRank;
    MPI_Comm_rank(comm, &myRank);
    std::cout << "Rank: " << myRank << std::endl;
    for(int i = 0; i < localSize; i++)
      std::cout << begin[i] << "\t";
    std::cout << std::endl;
    */

    //Sort 'cache' with size-changeable sorting function parallel_sort_changeable.
    int newLocalSize = parallel_sort_changeable(cache,localSize,comm);

    /*
    std::cout << "OK" ;
    */

    //Compute the configure for data transformation
    int *sendcnts, *sdispls, *recvcnts, *rdispls;
    sendcnts = new int[numProc];
    sdispls = new int[numProc];
    recvcnts = new int[numProc];
    rdispls = new int[numProc];
    Assign_Data(sendcnts, sdispls, recvcnts, rdispls, newLocalSize, localSize, comm);

    //Use alltoall to reallocate to [begin, end);
    MPI_Alltoallv(cache, sendcnts, sdispls, MPI_INTEGER, begin, recvcnts, rdispls, MPI_INTEGER, comm);

    delete [] cache;
    delete [] sendcnts;
    delete [] sdispls;
    delete [] recvcnts;
    delete [] rdispls;

    return ;
}


int parallel_sort_changeable(int* &begin, int arraysize, MPI_Comm comm)
{
    int myRank;
    int numProc;
    int totalArraySize;
    int localArraySize = arraysize;

    //Compute the size of whole array
    MPI_Allreduce(&localArraySize, &totalArraySize, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Comm_rank(comm, &myRank);
    MPI_Comm_size(comm, &numProc);

    /*
    int myRankInWorld;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankInWorld);
    std::cout << "Rank: " << myRankInWorld << std::endl;
    std::cout << "TotalSize:" << totalArraySize << ".\nLocalArray Get:\t";
    for(int i = 0; i < localArraySize; i++)
      std::cout << begin[i] << "\t";
    std::cout << std::endl;
    */


    /*************************************/
    /**** Processor more than elements*****/
    /*************************************/
    /*
     * construct new comm with 0,1,..,n-1;
     * The rest (n,n+1,...,p-1) rest
     */
    if (numProc > totalArraySize)
    {
        MPI_Comm newComm;
        if (myRank < totalArraySize)
        {
            MPI_Comm_split(comm, 0, myRank, &newComm);
            return parallel_sort_changeable(begin, arraysize, newComm);
        }
        else 
          return 0;
    }

    /*************************************/
    /******* One processor case **********/
    /*************************************/
    if (numProc == 1)
    {
        std::sort(begin, begin+arraysize);
        return arraysize;
    }

    /*************************************/
    /***** General Case with 1 < p <= n **/
    /*************************************/

    /** A: randomly choose the pivot between 0 and totalArraySize-1 **/

    //set the same seed (Numproc)

    //generate the pivot and broadcast it among comm
    int pivot_rank;
    int pivot;

    //Randomly generate the processor rank where the pivot locates.
    if (myRank == 0)
        pivot_rank = rand() % numProc;
    MPI_Bcast(&pivot_rank, 1, MPI_INTEGER, 0, comm);

    if (myRank == pivot_rank)
    {
        int pivot_index = rand() % block_decompose(totalArraySize, numProc, pivot_rank);
        pivot = *(begin+pivot_index);
    }
    MPI_Bcast(&pivot, 1, MPI_INTEGER, pivot_rank, comm);

    /** B: do the partition locally **/
    //the number of smaller-than-or-equal-to-pivot elements;
    int numLeft, numEqu;   
    Partition_Local(begin, begin + arraysize, pivot, numLeft, numEqu);
    int numTotalLeft, numTotalEqu;
    MPI_Allreduce(&numLeft, &numTotalLeft, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allreduce(&numEqu, &numTotalEqu, 1, MPI_INTEGER, MPI_SUM, comm);


    /** C: Deal the array with same entry case **/
    if (numTotalEqu == totalArraySize)
        return localArraySize;

    /*
    //SEQUENCE AFTER PIVOTING
    std::cout << "Pivot: " << pivot << ". LocalnumLeft:" << numLeft << ". GlobalnumLeft:" << numTotalLeft 
        << "\nArray After Pivoting:";
    for(int i = 0; i < localArraySize; i++)
      std::cout << begin[i] << "\t";
    std::cout << std::endl;
    */

    /** D: Assign new subproblems, calculate the destination of sending and receiving **/
    int *sendcnts, *sdispls, *recvcnts, *rdispls;
    sendcnts = new int[numProc];
    sdispls = new int[numProc];
    recvcnts = new int[numProc];
    rdispls = new int[numProc];

    int *tmp_sendcnts, *tmp_sdispls, *tmp_recvcnts, *tmp_rdispls;
    tmp_sendcnts = new int[numProc];
    tmp_sdispls = new int[numProc];
    tmp_recvcnts = new int[numProc];
    tmp_rdispls = new int[numProc];

    //Partition the numProc processor
    int numProcLeft,numProcRight;
    Partition_Proc(numProcLeft, numProcRight, numProc, numLeft, localArraySize - numLeft, comm);


    //New assigned local array size
    int newLocalSize;
    if (myRank < numProcLeft)
      newLocalSize = block_decompose(numTotalLeft, numProcLeft, myRank);
    else
      newLocalSize =block_decompose(totalArraySize - numTotalLeft, numProcRight, myRank - numProcLeft);

    /*
    //Partition of processor and new localsize output
    std::cout << "Partition of Processor-Left: " << numProcLeft << "\t Right: " << numProcRight << std::endl;
    std::cout << "NewLocalSize: " << newLocalSize << std::endl;
    */


    //Compute the configuration for the data transform
    if(myRank < numProcLeft)
    {
        //For the smaller-than-equal-to elements
        Assign_Data(tmp_sendcnts, tmp_sdispls, tmp_recvcnts, tmp_rdispls, numLeft, newLocalSize, comm);
        //For the larger-than-pivot elements
        Assign_Data(sendcnts, sdispls, recvcnts, rdispls, localArraySize - numLeft, 0, comm);
    }
    else
    {
        //For the smaller-than-equal-to elements
        Assign_Data(tmp_sendcnts, tmp_sdispls, tmp_recvcnts, tmp_rdispls, numLeft, 0, comm);
        //For the larger-than-pivot elements
        Assign_Data(sendcnts, sdispls, recvcnts, rdispls, localArraySize - numLeft, newLocalSize, comm);
    }

    //Combine the two configuration
    if(myRank < numProcLeft)
    {
        for(int i = 0; i < numProcLeft; i++)
        {
            if(sendcnts[i] != 0)
              std::cout << "EERROR IN CONFIGURATION OF ALLTOALL\n";
            sendcnts[i] = tmp_sendcnts[i];
        }
        for(int i = 0; i < numProc; i++)
          recvcnts[i] = tmp_recvcnts[i];
    }
    else
    {
        for(int i = 0; i < numProcLeft; i++)
          sendcnts[i] = tmp_sendcnts[i];
    }

    delete [] tmp_sendcnts;
    delete [] tmp_sdispls ;
    delete [] tmp_recvcnts;
    delete [] tmp_rdispls ;

    //Recalculate the displacement
    sdispls[0] = 0; rdispls[0] = 0;
    for(int i = 1; i < numProc; i++)
    {
        sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        rdispls[i] = rdispls[i-1] + recvcnts[i-1];
    }

    //Exchange Data based on the computed configuration 
    int* tmp = new int[newLocalSize];
    MPI_Alltoallv(begin, sendcnts, sdispls, MPI_INTEGER, tmp, recvcnts, rdispls, MPI_INTEGER, comm);

    delete [] begin;
    begin = new int[newLocalSize];
    for(int i = 0; i < newLocalSize; i++)
      begin[i] = *(tmp+i);
    delete [] tmp;

    /** F: Construct new communicator and sort **/
    //New comm
    MPI_Comm newComm;
    MPI_Comm_split(comm, (myRank < numProcLeft), myRank, &newComm);


    //Size of new comm
    int finalSize = parallel_sort_changeable(begin, newLocalSize, newComm);

    delete [] sendcnts;
    delete [] sdispls;
    delete [] recvcnts;
    delete [] rdispls;
    MPI_Comm_free(&newComm);

    return finalSize;
}


void Partition_Proc(int &procLeft, int &procRight, int numProc, int numLeft, int numRight, MPI_Comm comm)
{
    int TotalLeft, TotalRight;
    MPI_Allreduce(&numLeft, &TotalLeft, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allreduce(&numRight, &TotalRight, 1, MPI_INTEGER, MPI_SUM, comm);

    if(TotalLeft < TotalRight)
      procLeft = (int)floor(double(numProc * TotalLeft)/(TotalLeft + TotalRight));
    else
      procLeft = (int)ceil(double(numProc * TotalLeft)/(TotalLeft + TotalRight));
    procRight = numProc - procLeft;

    if(procLeft == 0 && TotalLeft != 0)
    {
        procLeft++;
        procRight--;
    }

    if(procRight == 0 && TotalRight != 0)
    {
        procRight++;
        procLeft--;
    }

    if((procRight == 0 && TotalRight != 0) || (procLeft == 0 && TotalLeft != 0))
      std::cout << "ERROR HAPPENED IN PARTITION OF PROCESSOR\n";

    return ;
}


void Partition_Local(int* begin, int* end, int pivot, int &numLeft, int& numEqn)
{
    int index = 0;
    int swap_tmp;
    numEqn = 0;
    for(int i = 0; i < end - begin; i++)
    {
        if(begin[i] <= pivot)
        {
            swap_tmp = begin[index];
            begin[index] = begin[i];
            begin[i] = swap_tmp;
            index++;
            if(begin[i] == pivot)
              numEqn++;
        }
    }
    numLeft = index;
    return ;
}

void Assign_Data(int* sendcnts, int* sdispls, 
            int* recvcnts, int* rdispls, int curSize, int objSize, MPI_Comm &comm)
{
    int myRank;
    int numProc;
    MPI_Comm_size(comm, &numProc);
    MPI_Comm_rank(comm, &myRank);

    int* prefix_curSize = new int[numProc];
    int* prefix_objSize = new int[numProc];
    int tmp;

    //Collect the prefix sum of present size and allgather them
    MPI_Scan(&curSize, &tmp, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allgather(&tmp, 1, MPI_INTEGER, prefix_curSize, 1, MPI_INTEGER, comm);

    //Collect the prefix sum of objective size and allgather them
    MPI_Scan(&objSize, &tmp, 1, MPI_INTEGER, MPI_SUM, comm);
    MPI_Allgather(&tmp, 1, MPI_INTEGER, prefix_objSize, 1, MPI_INTEGER, comm);

    int head, tail;

    //Find the index of processors [head, tail] to RECEIVE data from
    head = 0; tail = 0;
    if (myRank == 0)
      head = 0;
    else
      while(prefix_curSize[head] < prefix_objSize[myRank-1]) head++;

    tail = head;
    while(prefix_curSize[tail] < prefix_objSize[myRank]) tail++;

    //Setup the RECEIVE counts
    for(int i = 0; i < head; i++)
      recvcnts[i] = 0;

    if(tail == head)
      recvcnts[head] = objSize;
    else
    {
        //For head processor
        recvcnts[head] = prefix_curSize[head] - ((myRank > 0) ? prefix_objSize[myRank-1] : 0);
        //For processor between (head, tail)
        for(int i = head + 1; i < tail; i++)
          recvcnts[i] = prefix_curSize[i] - prefix_curSize[i-1];
        //For tail processor
        recvcnts[tail] = prefix_objSize[myRank] - prefix_curSize[tail - 1];
    }
    for(int i = tail + 1; i < numProc; i++) 
      recvcnts[i] = 0;

    //Find the index of processors [head, tail] to SEND data
    head = 0; tail = 0;
    if (myRank == 0)
      head = 0;
    else
      while(prefix_objSize[head] < prefix_curSize[myRank-1]) head++;

    tail = head;
    while(prefix_objSize[tail] < prefix_curSize[myRank]) tail++;

    //Setup the send counts
    for(int i = 0; i < head; i++)
      sendcnts[i] = 0;

    if(tail == head)
      sendcnts[head] = curSize;
    else
    {
        //For head processor
        sendcnts[head] = prefix_objSize[head] - ((myRank > 0) ? prefix_curSize[myRank-1] : 0);
        //For processor between (head, tail)
        for(int i = head + 1; i < tail; i++)
          sendcnts[i] = prefix_objSize[i] - prefix_objSize[i-1];
        //For tail processor
        sendcnts[tail] = prefix_curSize[myRank] - prefix_objSize[tail - 1];
    }
    for(int i = tail + 1; i < numProc; i++)
      sendcnts[i] = 0; 

    //SET UP THE DISPLACEMENT;
    sdispls[0] = 0; rdispls[0] = 0;
    for(int i = 1; i < numProc; i++)
    {
        sdispls[i] = sdispls[i-1] + sendcnts[i-1];
        rdispls[i] = rdispls[i-1] + recvcnts[i-1];
    }

    delete [] prefix_curSize;
    delete [] prefix_objSize;

    return ;
}


