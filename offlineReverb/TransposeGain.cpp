//
//  TransposeGain.c
//  MatrixGainTransposition
//
//  Created by Hans on 15/3/17.
//  Copyright © 2017 Hans. All rights reserved.
//

#include "TransposeGain.hpp"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <assert.h>



/*
 * X,Y, and Z are square matricies
 * 
 * calculates Z = X.Y
 */
void squareMatrixMultiply(float** X, float** Y, float** Z, size_t N){
    // initialize the output to zero
    for(size_t i=0; i<N; i++)
        memset(Z[i],0,sizeof(float)*N);
    
    // iterate on output row
    for(size_t i = 0; i<N; i++){
        // iterate on output column
        for(size_t j=0; j<N; j++){
            // iterate down input vectors
            for(size_t k=0; k<N; k++){
                Z[i][j] += X[i][k]*Y[k][j];
            }
        }
    }
}




/*
 * computes y = M.x
 *
 * @param M is a square matrix of size length
 * @param x is an array of length
 * @param y != x is an array of length
 *
 */
void matrixVectorMultiply(float** M, float* x, float* y, size_t length){
    // init the output to zero
    memset(y,0,sizeof(float)*length);
    
    for(size_t i=0; i<length; i++)
        for(size_t j=0; j<length; j++)
            y[i] += M[i][j] * x[j];
}





// hadamard transform for vectors of length 2
// does not work in place
void hadamard2(float* input, float* output){
    output[0] = input[0] + input[1];
    output[1] = input[0] - input[1];
}






// does not work in place
void hadamardTransform(float* input, float* output, int length){
    int partitionSize = length;
    
    float* temp1 = new float[length];
    float* temp2 = new float[length];
    
    float* tmpIn = input;
    float* tmpOut = temp1;
    bool tmpState = true;
    
    // iteratively calculated recursion
    while (partitionSize > 0) {
        int halfPartitionSize = partitionSize >> 1;
        for (int i = 0; i < length; i += partitionSize){
            // copy all the lower terms into place
            memcpy(tmpOut+i,tmpIn+i,sizeof(float)*halfPartitionSize);
            memcpy(tmpOut+i+halfPartitionSize,tmpIn+i,sizeof(float)*halfPartitionSize);
            // sum all the higher terms into place
            for (int j=i; j<halfPartitionSize+i; j++) {
                int idx2 = j+halfPartitionSize;
                tmpOut[j] += tmpIn[idx2];
                tmpOut[idx2] -= tmpIn[idx2];
            }
        }
        
        // swap temp buffers to avoid using the same memory for reading and writing.
        tmpIn = tmpOut;
        tmpState = !tmpState;
        if (tmpState) tmpOut = temp1;
        else tmpOut = temp2;
        partitionSize >>= 1;
    }
    
    
    // base case (size 2 transform)
    for (int i = 0; i < length; i += 2)
        hadamard2(tmpIn + i, output + i);
    
    // free the temp arrays
    free(temp1);
    free(temp2);
}



// allocate memory for a matrix of size length
float** mallocMatrix(size_t length){
    float** M;
    
    // malloc first dimension
    M = new float*[length];
    
    // malloc second dimension
    for(size_t i=0; i<length; i++)
        M[i] = new float[length];
    
    return M;
}






// free memory for matrix of size length
void freeMatrix(float** M, size_t length){
    // free second dimension
    for(size_t i=0; i<length; i++)
        free(M[i]);
    
    // free first dimension
    free(M);
}






// make a diagonal matrix of vector v and store it
// into M, of size length.
void toDiagonalMatrix(float* v, float** M, size_t length){
    
    // initialize the output to zero
    for(size_t i=0; i<length; i++)
        memset(M[i],0,sizeof(float)*length);
    
    // put V along the diagonal
    for(size_t i=0; i<length; i++)
        M[i][i] = v[i];
}






void identityMatrix(float**M, size_t length){
    // initialize the output to zero
    for(size_t i=0; i<length; i++)
        memset(M[i],0,sizeof(float)*length);
    
    // put 1 along the diagonal
    for(size_t i=0; i<length; i++)
        M[i][i] = 1.0f;
}






// N = s*M
void matrixScalarMul(float** M, float s, float** N, size_t length){
    for(size_t i=0; i<length; i++)
        for(size_t j=0; j<length; j++)
            N[i][j]=M[i][j]*s;
}






void makeHadamardMatrix(float**M, size_t length){
    float e = 1.0f/sqrt((float)length);
    
    // create an identity matrix
    identityMatrix(M,length);
    
    // scale the identity matrix down (because our hadamard transform
    // requires that to make it unitary)
    matrixScalarMul(M,e,M,length);
    
    // create a temp array
    float* t = new float[length];
    
    // take the hadamard transform of each row of the scaled identity Matrix
    for(size_t i=0; i<length; i++){
        // transform into temp array
        hadamardTransform(M[i],t,(int)length);
        
        // copy back to M
        memcpy(M[i], t, sizeof(float)*length);
    }
    
    // free the temp array
    free(t);
}







void printMatrix(float** M, size_t length){
    printf("{");
    
    for(size_t i=0; i<length; i++){
        printf("{");
        
        for(size_t j=0; j<length; j++){
            if(j==length-1)
                printf("%f", M[i][j]);
            else
                printf("%f, ", M[i][j]);
        }
        
        if(i==(length-1))
            printf("}\n");
        else
            printf("},\n");
    }
    printf("}\n");
}



bool isPowerOfTwo (size_t x)
{
    return ((x != 0) && ((x & (~x + 1)) == x));
}




// returns the L2 norm of x
float norm(float* x, size_t length){
    float sum = 0.0f;
    
    for(size_t i=0; i<length; i++)
        sum += x[i]*x[i];
    
    return sqrt(sum);
}



/*
 * returns (energy with upsilon at input) / (energy with upsilon at output)
 *
 *      Norm[M.U.ones] / Norm[U.M.ones]
 *
 * to get the gain coefficient to compensate for this, take the square root
 * of the output.
 *
 * @param upsilon - an array representing the output tap coefficients
 *                  for an imagingary FDN 
 * @length        - the length of upsilon
 */
float transposeEnergyRatio(float* upsilon, size_t length){
    assert(isPowerOfTwo(length));
    
    float** M; // Mixing Matrix (Hadamard)
    float** U; // Upsilon array as diagonal matrix
    float** Z; // output of matrix multiplication
    
    M = mallocMatrix(length);
    U = mallocMatrix(length);
    Z = mallocMatrix(length);
    
    makeHadamardMatrix(M, length);
    
    // print M (should be the identity matrix)
    //printMatrix(M, length);
    
    // confirm that M.M is the identity matrix (i.e. M is unitary)
    //squareMatrixMultiply(M, M, Z, length);
    
    // print M.M (should be the identity matrix)
    //printMatrix(Z, length);
    
    // put the vector upsilon on the diagonal of U
    toDiagonalMatrix(upsilon, U, length);
    
    // create a vector of ones to represent equal energy input on all channels
    float* ones = new float[length];
    for(size_t i=0;i<length; i++) ones[i]=1.0f;
    
    // create an empty vector for output
    float* y = new float[length];
    
    // Z = M.U
    squareMatrixMultiply(M, U, Z, length);
    
    // y = Z.ones
    matrixVectorMultiply(Z, ones, y, length);
    
    float afterNorm = norm(y, length);
    
    // redo the matrix multiplication with upsilon on the left side
    // Z = U.M
    squareMatrixMultiply(U, M, Z, length);
    
    // y = Z.ones
    matrixVectorMultiply(Z, ones, y, length);
    
    float beforeNorm = norm(y, length);
    
    
    // find the ration of energy before / after
    float energyChange = (afterNorm*afterNorm)/ (beforeNorm*beforeNorm) ;
    
    // free memory
    free(ones);
    free(y);
    freeMatrix(M, length);
    freeMatrix(U, length);
    freeMatrix(Z, length);
    
    return energyChange;
}
