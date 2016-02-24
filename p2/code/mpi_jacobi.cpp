/**
 * @file    mpi_jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements MPI functions for distributing vectors and matrixes,
 *          parallel distributed matrix-vector multiplication and Jacobi's
 *          method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */

#include "mpi_jacobi.h"
#include "jacobi.h"
#include "utils.h"

#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <vector>

void distribute_vector(const int n, double* input_vector, double** local_vector, MPI_Comm comm)
{
    //Information of present processor.
    int MyRank;
    int MyCoord[2];
    MPI_Comm_rank(comm, &MyRank);
    MPI_Cart_coords(comm, MyRank, 2, MyCoord);

    //1.Create column communicator to do the scatter operation only for first column
    MPI_Comm comm_col;
    int remain_dims[2]={true,false};
    MPI_Cart_sub(comm,remain_dims,&comm_col);

    MPI_Group CartGroup, ColGroup;
    MPI_Comm_group(comm, &CartGroup);
    MPI_Comm_group(comm_col, &ColGroup); 

    //If not in the first column, no work needed.
    if (MyCoord[1] != 0)
    {
        MPI_Group_free(&CartGroup);
        MPI_Group_free(&ColGroup);
        MPI_Comm_free(&comm_col);
        return ;
    }

    //2.Find the root rank at (0,0);
    int RootRank;
    int RootCoord[] = {0,0};
    MPI_Cart_rank(comm, RootCoord, &RootRank);

    //3.Set the sendcounts and displs
    int dims[2];
    int period[2];
    int coords[2];
    MPI_Cart_get(comm,2,dims,period,coords);

    int* sendcounts = new int[dims[0]];
    int* displs = new int[dims[0]];
    displs[0] = 0;
    sendcounts[0] = block_decompose(n, dims[0], 0);

    for(int i = 1; i < dims[0]; i++)
    {    
        sendcounts[i] = block_decompose(n, dims[0], i);
        displs[i] = displs[i-1] + sendcounts[i-1]; 
    }

    //4.MPI_SCATTERV.
    int localsize = block_decompose(n, dims[0], MyCoord[0]);
    (*local_vector) = new double[localsize];
    int Comm_ColRootRank;
    MPI_Group_translate_ranks(CartGroup, 1, &RootRank, ColGroup, &Comm_ColRootRank);
    MPI_Scatterv(input_vector, sendcounts, displs, MPI_DOUBLE, *local_vector, localsize, MPI_DOUBLE, Comm_ColRootRank, comm_col);

    delete sendcounts;
    delete displs;
    MPI_Group_free(&CartGroup);
    MPI_Group_free(&ColGroup);
    MPI_Comm_free(&comm_col);
    return ;
}


// gather the local vector distributed among (i,0) to the processor (0,0)
void gather_vector(const int n, double* local_vector, double* output_vector, MPI_Comm comm)
{
    //Information of present processor.
    int MyRank;
    int MyCoord[2];
    MPI_Comm_rank(comm, &MyRank);
    MPI_Cart_coords(comm, MyRank, 2, MyCoord);

    //1.Create column communicator to do the scatter operation only for first column
    MPI_Comm comm_col;
    int remain_dims[2]={true,false};
    MPI_Cart_sub(comm,remain_dims,&comm_col);

    MPI_Group CartGroup, ColGroup;
    MPI_Comm_group(comm, &CartGroup);
    MPI_Comm_group(comm_col, &ColGroup); 

    //If not in the first column, no work needed.
    if (MyCoord[1] != 0)
    {
        MPI_Group_free(&CartGroup);
        MPI_Group_free(&ColGroup);
        MPI_Comm_free(&comm_col);
        return ;
    }

    //2.Find the root rank at (0,0);
    int RootRank;
    int RootCoord[] = {0,0};
    MPI_Cart_rank(comm, RootCoord, &RootRank);

    // Information of the cartesian communicator
    int dims[2];
    int period[2];
    int coords[2];
    MPI_Cart_get(comm,2,dims,period,coords);

    // Build the recvcount and displs array for mpi_gatherv
    int* recvcount = NULL;
    int* displs = NULL;
    if(MyRank == RootRank)
    {
        recvcount = new int[dims[0]];
        displs = new int[dims[0]];
        recvcount[0] = block_decompose(n, dims[0], 0);
        displs[0] = 0;

        for(int i = 1; i < dims[0]; i++)
        {
            recvcount[i] = block_decompose(n, dims[0], i);
            displs[i] = displs[i-1] + recvcount[i-1];
        }
    }

    // MPI_GATHERV
    int localsize = block_decompose(n, dims[0], MyCoord[0]);
    int Comm_ColRootRank;
    MPI_Group_translate_ranks(CartGroup, 1, &RootRank, ColGroup, &Comm_ColRootRank);
    MPI_Gatherv(local_vector, localsize, MPI_DOUBLE, output_vector, recvcount, displs, MPI_DOUBLE, Comm_ColRootRank, comm_col);

    delete recvcount;
    delete displs;
    MPI_Group_free(&CartGroup);
    MPI_Group_free(&ColGroup);
    MPI_Comm_free(&comm_col);
    return ;
}

void distribute_matrix(const int n, double* input_matrix, double** local_matrix, MPI_Comm comm)
{
    //Step1. Scatter n/q rows of the matrix from (0,0) process to (i,0)
    //Step2. Use a for-loop to scatter the n/q columns in each row of the submatrix in (i,0) to (i,j);

    //Find the root rank at (0,0);
    int RootRank;
    int RootCoord[] = {0,0};
    MPI_Cart_rank(comm, RootCoord, &RootRank);

    //Find the basic information of the comm
    int dims[2];
    int period[2];
    int coords[2];
    MPI_Cart_get(comm,2,dims,period,coords);

    //Find the information of this processor.
    int MyRank;
    int MyCoord[2];
    MPI_Comm_rank(comm, &MyRank);
    MPI_Cart_coords(comm, MyRank, 2, MyCoord);

    //Located-Column Communicator
    MPI_Comm comm_col;
    int remain_dims[2]={true,false};
    MPI_Cart_sub(comm,remain_dims,&comm_col);

    //Located-Row Communicator
    MPI_Comm comm_row;
    remain_dims[0]= false; remain_dims[1]= true;
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    //Corresponding groups
    MPI_Group CartGroup, ColGroup, RowGroup;
    MPI_Comm_group(comm, &CartGroup);
    MPI_Comm_group(comm_col, &ColGroup); 
    MPI_Comm_group(comm_row, &RowGroup); 

    //Step 1.  in the first column, scatter the rows of the matrix.
    double* intermediate_matrix = NULL;
    if (MyCoord[1] == 0)
    {
        //Set the sendcounts and displs
        int* sendcounts = new int[dims[0]];
        int* displs = new int[dims[0]];

        displs[0] = 0;
        sendcounts[0] = n * block_decompose(n, dims[0], 0);

        for(int i = 1; i < dims[0]; i++)
        {    
            sendcounts[i] = n * block_decompose(n, dims[0], i);
            displs[i] = displs[i-1] + sendcounts[i-1]; 
        }

        //4.MPI_SCATTERV.
        int localsize = n * block_decompose(n, dims[0], MyCoord[0]);
        intermediate_matrix = new double[localsize];
        int Comm_RootRank;
        MPI_Group_translate_ranks(CartGroup, 1, &RootRank, ColGroup, &Comm_RootRank);
        MPI_Scatterv(input_matrix, sendcounts, displs, MPI_DOUBLE, intermediate_matrix, localsize, MPI_DOUBLE, Comm_RootRank, comm_col);

        delete sendcounts;
        delete displs;
    }
    MPI_Barrier(comm);

    //Step 2. Scatter the columns of the submatrix

    //The submatrix size in each processor
    int local_rows = block_decompose(n, dims[0], MyCoord[0]);
    int local_cols = block_decompose(n, dims[1], MyCoord[1]);
    (*local_matrix) = new double[local_rows * local_cols];

    //The row root rank at the same 
    int RowRootRank;
    int RowRootCoord[2] = {MyCoord[0], 0};
    MPI_Cart_rank(comm, RowRootCoord, &RowRootRank);

    //Set the sendcounts and displs for each row scatter
    int* sendcounts = new int[dims[1]];
    int* displs = new int[dims[1]];

    displs[0] = 0;
    sendcounts[0] = block_decompose(n, dims[1], 0);

    for(int i = 1; i < dims[0]; i++)
    {    
        sendcounts[i] = block_decompose(n, dims[1], i);
        displs[i] = displs[i-1] + sendcounts[i-1]; 
    }


    //Scatter each row in (i, 0) processor
    int Comm_RowRootRank;
    MPI_Group_translate_ranks(CartGroup, 1, &RowRootRank, RowGroup, &Comm_RowRootRank);
    for(int i = 0; i < local_rows; i++)
      MPI_Scatterv((intermediate_matrix + i*n), sendcounts, displs, MPI_DOUBLE, (*local_matrix + i*local_cols), local_cols, MPI_DOUBLE, Comm_RowRootRank, comm_row);

    delete sendcounts;
    delete displs;
    delete intermediate_matrix;
    MPI_Comm_free(&comm_row);
    MPI_Comm_free(&comm_col);
    MPI_Group_free(&CartGroup);
    MPI_Group_free(&RowGroup);
    MPI_Group_free(&ColGroup);
    return ;
}


void transpose_bcast_vector(const int n, double* col_vector, double* row_vector, MPI_Comm comm)
{
    //Require the col_vector, row_vector have pre-assigned space. 
    //Step 1: Send the vector in (i,0) to (i,i);
    //Step 2: Broadcast the vector in (i,i) to the ith column.

    //Information of the comm
    int dims[2];
    int period[2];
    int coords[2];
    MPI_Cart_get(comm,2,dims,period,coords);

    //Information of this processor.
    int MyRank;
    int MyCoord[2];
    MPI_Comm_rank(comm, &MyRank);
    MPI_Cart_coords(comm, MyRank, 2, MyCoord);

    //Information of the first-column processor(with the same row)
    int RowRootRank;
    int RowRootCoord[2] = {MyCoord[0], 0};
    MPI_Cart_rank(comm, RowRootCoord, &RowRootRank);

    //Information of the diagonal processor(with the same column & row)
    int ColDiagRootRank;
    int ColDiagRootCoord[2] = {MyCoord[1], MyCoord[1]};
    MPI_Cart_rank(comm, ColDiagRootCoord, &ColDiagRootRank);

    int RowDiagRootRank;
    int RowDiagRootCoord[2] = {MyCoord[0], MyCoord[0]};
    MPI_Cart_rank(comm, RowDiagRootCoord, &RowDiagRootRank);


    //Information of the vector;
    int col_size = block_decompose(n, dims[0], MyCoord[0]);
    int row_size = block_decompose(n, dims[0], MyCoord[1]);


    //Step 1: Send the vector in (i, 0) to (i,i)
    if (MyCoord[0] == 0 && MyCoord[1] == 0)
    {
        for(int i = 0; i < row_size; i++)
          row_vector[i] = col_vector[i];
    }
    else{
        if (MyCoord[1] == 0)
          MPI_Send(col_vector, col_size, MPI_DOUBLE, RowDiagRootRank, 1, comm);
        if (MyCoord[0] == MyCoord[1])
          MPI_Recv(row_vector, row_size, MPI_DOUBLE, RowRootRank, 1, comm, MPI_STATUS_IGNORE);
    }

    //Step 2: Broadcast over the column
    MPI_Comm comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    MPI_Group CartGroup, ColGroup;
    MPI_Comm_group(comm, &CartGroup); 
    MPI_Comm_group(comm_col, &ColGroup); 

    int Comm_DiagRank;
    MPI_Group_translate_ranks(CartGroup, 1, &ColDiagRootRank, ColGroup, &Comm_DiagRank);

    MPI_Barrier(comm);
    MPI_Bcast(row_vector, row_size, MPI_DOUBLE, Comm_DiagRank, comm_col);

    MPI_Comm_free(&comm_col);
    MPI_Group_free(&CartGroup);
    MPI_Group_free(&ColGroup);
    return ;
}


void distributed_matrix_vector_mult(const int n, double* local_A, double* local_x, double* local_y, MPI_Comm comm)
{
    //Step 1: First trans_bcast_vector of the given local_x at the first column
    //Step 2: Do the submatrix * subvector in each processor
    //Step 3: MPI_Reduce to get the sum over each row and store the result at the local_y at the first column

    //Information of the comm
    int dims[2];
    int period[2];
    int coords[2];
    MPI_Cart_get(comm,2,dims,period,coords);

    //Information of this processor.
    int MyRank;
    int MyCoord[2];
    MPI_Comm_rank(comm, &MyRank);
    MPI_Cart_coords(comm, MyRank, 2, MyCoord);

    //Information of the first-column processor(with the same row)
    int RowRootRank;
    int RowRootCoord[2] = {MyCoord[0], 0};
    MPI_Cart_rank(comm, RowRootCoord, &RowRootRank);

    //Step 1;
    double* subvec;
    int vecsize = block_decompose(n, dims[0], MyCoord[1]);
    subvec = new double[vecsize];
    transpose_bcast_vector(n, local_x, subvec, comm);

    //Step 2:
    double* result_y; //The result in each processor. 
    int nrows = block_decompose(n, dims[0], MyCoord[0]);
    int ncols = block_decompose(n, dims[1], MyCoord[1]);
    result_y = new double[nrows];


    MPI_Barrier(comm);
    for(int i = 0; i < nrows; i++)
    {
        result_y[i] = 0;
        for(int j = 0; j < ncols; j++)
          result_y[i] += local_A[i*ncols + j] * subvec[j];
    }

    //Step 3:
    MPI_Comm comm_row;
    int remain_dims[2] = {false, true};
    MPI_Cart_sub(comm, remain_dims, &comm_row);

    MPI_Group CartGroup, RowGroup;
    MPI_Comm_group(comm, &CartGroup);
    MPI_Comm_group(comm_row, &RowGroup); 

    int Comm_RowRootRank;
    MPI_Group_translate_ranks(CartGroup, 1, &RowRootRank, RowGroup, &Comm_RowRootRank);

    for(int i = 0; i < nrows; i++)
      MPI_Reduce(result_y + i, local_y + i, 1, MPI_DOUBLE, MPI_SUM, Comm_RowRootRank, comm_row);

    MPI_Group_free(&CartGroup);
    MPI_Group_free(&RowGroup);
    MPI_Comm_free(&comm_row);
    delete subvec;
    delete result_y;
    return ;
}

// Solves Ax = b using the iterative jacobi method
void distributed_jacobi(const int n, double* local_A, double* local_b, double* local_x,
            MPI_Comm comm, int max_iter, double l2_termination)
{
    //Information of the comm
    int dims[2];
    int period[2];
    int coords[2];
    MPI_Cart_get(comm,2,dims,period,coords);

    //Information of this processor.
    int MyRank;
    int MyCoord[2];
    MPI_Comm_rank(comm, &MyRank);
    MPI_Cart_coords(comm, MyRank, 2, MyCoord);

    //Information of the first-column processor(with the same row)
    int RowRootRank;
    int RowRootCoord[2] = {MyCoord[0], 0};
    MPI_Cart_rank(comm, RowRootCoord, &RowRootRank);

    //Information of the (0,0) processor 
    int RootRank;
    int RootCoord[2] = {0, 0};
    MPI_Cart_rank(comm, RootCoord, &RootRank);

    //Information of the diagonal processor(with the same column)
    int DiagRootRank;
    int DiagRootCoord[2] = {MyCoord[0], MyCoord[0]};
    MPI_Cart_rank(comm, DiagRootCoord, &DiagRootRank);

    //Information of the submatrix and subvector
    int nrows = block_decompose(n, dims[0], MyCoord[0]);


    //Store the diagonal elements in the first column processor.
    double* diag = NULL;
    if(MyCoord[1] == 0)
    {
        diag = new double[nrows];
        if (MyRank != DiagRootRank)
        {
            for(int i = 0; i < nrows; i++)
              MPI_Recv(diag+i, 1, MPI_DOUBLE, DiagRootRank, MPI_ANY_TAG, comm, MPI_STATUS_IGNORE);  
        }
        else
        {
            for(int i = 0; i < nrows; i++)
              diag[i] = local_A[i*nrows + i];
        }
        //initilize the value.
        for (int i = 0; i < nrows; i++)
          local_x[i] = 0;
    }
    else if (MyRank == DiagRootRank)
    {
        for(int i = 0; i < nrows; i++)
          MPI_Send((local_A + i*nrows + i), 1, MPI_DOUBLE, RowRootRank, 1, comm);
    }

    //Store the result of A * x in the first column processors. 
    double* Ax = NULL;  
    if (MyCoord[1] == 0)
      Ax = new double[nrows];

    //Column Comunicator and its corresponding ranks.
    MPI_Comm comm_col;
    int remain_dims[2] = {true, false};
    MPI_Cart_sub(comm, remain_dims, &comm_col);

    MPI_Group CartGroup, ColGroup;
    MPI_Comm_group(comm, &CartGroup);
    MPI_Comm_group(comm_col, &ColGroup); 

    int Comm_RootRank;
    MPI_Group_translate_ranks(CartGroup, 1, &RootRank, ColGroup, &Comm_RootRank);

    //Jacobi Iteration.
    double local_res;
    double residual;
    int itercount = 0;

    while(1)
    {
        //Compute the residual norm
        distributed_matrix_vector_mult(n, local_A, local_x, Ax, comm);

        if (MyCoord[1] == 0)
        {
            local_res = 0;
            for(int i = 0; i < nrows; i++)
              local_res += (Ax[i] -local_b[i])*(Ax[i] - local_b[i]);
            MPI_Reduce(&local_res, &residual, 1, MPI_DOUBLE, MPI_SUM, Comm_RootRank, comm_col);

        }
        if (MyRank == RootRank)
            residual = sqrt(residual);

        MPI_Bcast(&residual, 1, MPI_DOUBLE, RootRank, comm);

        //Judge the termination condition and compute new iteration
        if (itercount < max_iter){
            if(residual < l2_termination){
                std::cout<< "Iteration complete with given l2 norm constriant."<<std::endl;
                break;
            }//else not converge within the residual range
            else {
                //Compute the next iteration, still stored in local_x;
                if(MyCoord[1] == 0)
                {
                    for(int i = 0; i < nrows; i++)
                      local_x[i] = (- Ax[i] +  local_b[i]  + diag[i]*local_x[i])/diag[i]; 
                }
                itercount++;
            }
        }
        else if (residual < l2_termination)
        {
            std::cout<< "Iteration Complete with given l2 norm constriant."<<std::endl;
            break;
        }
        else //not converged
        {
            std::cout<< "Not converged after" << max_iter<<"iteration" << std::endl;
            break;
        }
    }

    delete Ax;
    delete diag;
    MPI_Comm_free(&comm_col);
    MPI_Group_free(&CartGroup);
    MPI_Group_free(&ColGroup);
    return ;
}


// wraps the distributed matrix vector multiplication
void mpi_matrix_vector_mult(const int n, double* A,
            double* x, double* y, MPI_Comm comm)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_x = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &x[0], &local_x, comm);
    //MPI_Barrier(comm);
    //std::cout << "ALLGOOD\n";

    // allocate local result space
    double* local_y = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_matrix_vector_mult(n, local_A, local_x, local_y, comm);

    // gather results back to rank 0
    gather_vector(n, local_y, y, comm);
}

// wraps the distributed jacobi function
void mpi_jacobi(const int n, double* A, double* b, double* x, MPI_Comm comm,
            int max_iter, double l2_termination)
{
    // distribute the array onto local processors!
    double* local_A = NULL;
    double* local_b = NULL;
    distribute_matrix(n, &A[0], &local_A, comm);
    distribute_vector(n, &b[0], &local_b, comm);

    // allocate local result space
    double* local_x = new double[block_decompose_by_dim(n, comm, 0)];
    distributed_jacobi(n, local_A, local_b, local_x, comm, max_iter, l2_termination);

    // gather results back to rank 0
    gather_vector(n, local_x, x, comm);
}
