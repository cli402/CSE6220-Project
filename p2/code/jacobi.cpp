/**
 * @file    jacobi.cpp
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Implements matrix vector multiplication and Jacobi's method.
 *
 * Copyright (c) 2014 Georgia Institute of Technology. All Rights Reserved.
 */
#include "jacobi.h"
#include "utils.h"
#include <stdexcept>
#include <cstdlib>

// my implementation:
#include <iostream>
#include <cmath>

// Calculates y = A*x for a square n-by-n matrix A, and n-dimensional vectors x
// and y
void matrix_vector_mult(const int n, const double* A, const double* x, double* y)
{
    for(int i=0;i<n;i++){
        y[i] = 0;
        for (int j=0;j<n;j++)
        {
            y[i]+=A[i*n + j] * x[j];
        }
    } 
}

// Calculates y = A*x for a n-by-m matrix A, a m-dimensional vector x
// and a n-dimensional vector y
void matrix_vector_mult(const int n, const int m, const double* A, const double* x, double* y)
{
    for (int i=0;i<n;i++){
        y[i] = 0;
        for (int j =0;j<m;j++){
            y[i]+=A[i*m + j]*x[j];
        }
    }
}

// implements the sequential jacobi method
void jacobi(const int n, double* A, double* b, double* x, int max_iter, double l2_termination)
{
    //Use two temporary array to finish the Jacobi iteration
    double** iter = new double*[2];
    iter[0] = new double[n];
    iter[1] = new double[n];

    //Set initial solution as x0 = 0;
    for (int i = 0; i < n; i++) iter[0][i] = 0; 

    //Set the count and residual norm for the iteration
    int count = 0;
    double norm;

    while (1){
        //Compute the residual norm. 
        int curr = count % 2;
        norm = 0;
        for(int i = 0; i < n; i++)
        {
            double bb = 0;
            for(int j = 0; j < n; j++)
              bb += A[i*n + j]*iter[curr][j];
            norm += (bb - b[i])*(bb - b[i]);
        }
        norm = sqrt(norm);

        //Decide whether to terminate the loop.
        if (count < max_iter){
            if(norm < l2_termination){
                std::cout<< "Iteration complete with given l2 norm constriant."<<std::endl;
                break;
            }//else not converge within the residual range
            else {
                //index of the current vector and next vector.
                int curr = count % 2;
                int next = (count+1) % 2;

                //Do one Jacobi iteration step.
                for(int i = 0; i < n; i++)
                {
                    iter[next][i] = b[i];
                    for(int j = 0; j < i; j++)
                      iter[next][i] -= A[i*n + j]*iter[curr][j];
                    for(int j = i+1; j < n; j++)
                      iter[next][i] -= A[i*n + j]*iter[curr][j];
                    iter[next][i] /= A[i*n+i];
                }
                count++;
            }
        }
        else if (norm < l2_termination)
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

    //copy the result to x.
    for(int i = 0; i < n; i++)
      x[i] = iter[count%2][i];
    delete iter[0];
    delete iter[1];
    delete iter;
    return ;
}
