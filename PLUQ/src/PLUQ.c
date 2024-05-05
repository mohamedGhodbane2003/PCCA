#include "PLUQ.h"

// Crout method with AVX2 instructions
void pluq_crout_avx2(Matrix* A, int** P, int** Q, int* rank, int p) {
    int m = A->rows;
    int n = A->cols;

    *P = createRange(m);
    *Q = createRange(n);

    int matrixRank = 0;
    int nullity = 0;
    // Allocate temporary arrays for AVX2 operations
    int *tmp1 = malloc((n - matrixRank) * sizeof(int));
    int *tmp2 = malloc(((m - nullity) - (matrixRank + 1)) * sizeof(int));

    while(matrixRank + nullity < m) {
        // Compute vector-matrix multiplication using AVX2 instructions
        vectorMatrixMultiplication_avx2(A->data, n, matrixRank, tmp1, p);
        // Update the matrix based on the computed results using AVX2 instructions
        updateMatrix1_avx2(A->data, matrixRank, tmp1, n, p);

        int pivot = matrixRank;
        while(pivot < n && A->data[matrixRank, pivot] == 0)
            pivot++; 
        if(pivot == n) {
            rowRotation(A, matrixRank, *P);
            nullity++;
        }else{
            int invpivot = inverse(A->data[matrixRank * n + pivot], p);
            // Compute scalar products for row elimination using AVX2 instructions
            for(int i = matrixRank + 1; i < m - nullity; i++){
                tmp2[i - matrixRank - 1] = scalar_product_avx2(A->data, matrixRank,i, pivot, n, p);
            }
            
            updateMatrix2_avx2(A->data, matrixRank, tmp2, n, m, nullity, pivot, invpivot, p);
            colRotation(A, matrixRank, pivot+1, *Q);
            matrixRank++;
        }
    }
    *rank = matrixRank;
}

// Perform the PLUQ decomposition using the Crout method
void pluq_crout(Matrix* A, int** P, int** Q, int* rank, int p) {
    int m = A->rows;
    int n = A->cols;

    *P = createRange(m);
    *Q = createRange(n);

    int matrixRank = 0;
    int nullity = 0;

    while(matrixRank + nullity < m) {
        for (int j = matrixRank; j < n; j++) {// Perform the scalar product
                int sum = 0;
                for (int k = 0; k < matrixRank; k++) {
                    sum = add(sum, mult(A->data[matrixRank * n + k], A->data[k * n + j], p), p);
                }
                 // Subtract the sum from the current element
                A->data[matrixRank * n + j] = sub(A->data[matrixRank * n + j], sum, p);
            }
        int pivot = matrixRank;// Initialize pivot index
        // Find the pivot element in the current row
        while(pivot < n && A->data[matrixRank, pivot] == 0)
            pivot++;
        
        if(pivot == n) { // If pivot column is beyond the matrix boundary, perform row rotation and increment nullity
            rowRotation(A, matrixRank, *P); // 'matrixRank' to the last row (m-1) and update P
            nullity++;
        }else{
            int invpivot = inverse(A->data[matrixRank * n + pivot], p);
            for(int i = matrixRank + 1; i < m - nullity; i++){// Perform the scalar product
                int tmp = 0;
                for(int j = 0; j < matrixRank; j++){
                    tmp = add(tmp, mult(A->data[i * n + j], A->data[j * n + pivot], p), p);
                }
                // Subtract the sum from the current element
                A->data[i * n + pivot] = sub(A->data[i * n + pivot], tmp, p);
                // Scale the pivot column
                A->data[i * n + pivot] = mult(A->data[i * n + pivot], invpivot, p);
            }
            colRotation(A, matrixRank, pivot+1, *Q); // Perform column rotation 'matrixRank' to 'pivot+1'
            matrixRank++;
        }
    }
    *rank = matrixRank;
}

// Perform the PLUQ decomposition using AVX2 instructions
void pluq_inplace_avx2(Matrix* A, int** P, int** Q, int* rank, int p) {
    int m = A->rows;
    int n = A->cols;

    *P = createRange(m);
    *Q = createRange(n);

    int matrixRank = 0;
    int nullity = 0;

    while (matrixRank + nullity < m) {
        int pivot = matrixRank;
        while (pivot < n && A->data[matrixRank * n + pivot] == 0)
            pivot += 1;
        if (pivot == n) {
            rowRotation(A, matrixRank, *P);
            nullity++;
        } else {
            int inv = inverse(A->data[matrixRank * n + matrixRank], p);
            if (pivot != matrixRank)
                colTransposition(A, matrixRank, pivot, *Q);
            // Initialize AVX2 vector vu with the reciprocal of the parameter p as double precision floating-point values
            __m256d vu = _mm256_set1_pd(1.0 / p);
            // Initialize AVX2 vector vp with the parameter p as double precision floating-point values
            __m256d vp = _mm256_set1_pd(p);
            // Initialize AVX2 vector vp_128 with the parameter p as 32-bit integer values
            __m128i vp_128 = _mm_set1_epi32(p);
            for (int k = matrixRank + 1; k < m; k++) {
                A->data[k * n + matrixRank] = mult(A->data[k * n + matrixRank], inv, p);
                // Perform row elimination using AVX2 instructions
                rows_elimination_avx2(A->data, n, matrixRank,A->data[k * n + matrixRank], p ,vp ,vu, vp_128 ,k);
            }
            matrixRank++;
        }
    }
    *rank = matrixRank;
}

// Perform the PLUQ decomposition of a matrix in place
void pluq_inplace(Matrix* A, int** P, int** Q, int* rank, int p) {
    int m = A->rows;
    int n = A->cols;

    *P = createRange(m);
    *Q = createRange(n);

    int matrixRank = 0;
    int nullity = 0;

    while (matrixRank + nullity < m) {
        int pivot = matrixRank; // Initialize the pivot index for the current iteration
        // Find the pivot column by searching for the first nonzero element in the pivot row
        while (pivot < n && A->data[matrixRank * n + pivot] == 0)
            pivot += 1;
        if (pivot == n) {// If the pivot column is beyond the matrix boundary, perform row rotation and increment nullity
            rowRotation(A, matrixRank, *P); // 'matrixRank' to the last row (m-1) and update P
            nullity++;
        } else {
            // Compute the inverse of the pivot element
            int inv = inverse(A->data[matrixRank * n + matrixRank], p);
            // If the pivot column is not the same as the pivot row, perform column transposition and update Q
            if (pivot != matrixRank)
                colTransposition(A, matrixRank, pivot, *Q);
            // Eliminate elements below the pivot
            for (int k = matrixRank + 1; k < m; k++) {
                // Compute the multiplier for row reduction
                A->data[k * n + matrixRank] = mult(A->data[k * n + matrixRank], inv, p);
                // Update elements in the current row using row reduction
                for (int j = matrixRank + 1; j < n; j++)
                    A->data[k * n + j] = sub(A->data[k * n + j], mult(A->data[k * n + matrixRank], A->data[matrixRank * n + j], p), p);
            }
            matrixRank++;
        }
    }
    *rank = matrixRank;
}
// Perform PLUQ decomposition on the input matrix A and generate permutation matrices P and Q and the decomposed matrix LU
void PLUQ(Matrix A, int** P, Matrix* LU, int** Q, int* rank, int p, int algo) {
    *LU = createMatrix(A.rows, A.cols);// Create a matrix LU with the same dimensions as A
    copyMatrix(A, LU); // Copy the content of matrix A to LU
    switch(algo){
        case 1: pluq_inplace(LU, P, Q, rank, p);
        break;
        case 2: pluq_inplace_avx2(LU, P, Q, rank, p);
        break;
        case 3: pluq_crout(LU, P, Q, rank, p);
        break;
        case 4: pluq_crout_avx2(LU, P, Q, rank, p);
        break;
        default : break;
    }
}

// Expand the decomposed matrix LU into lower triangular matrix L and upper triangular matrix U.
void expand_PLUQ(Matrix LU, int rank, Matrix* L, Matrix* U) {
    int m = LU.rows;
    int n = LU.cols;

    *L = zerosMatrix(m, m); 
    *U = zerosMatrix(m, n); 

    for (int j = 0; j < rank; j++) {// Populate the lower triangular part of L and diagonal elements with appropriate values from LU
        for (int i = j + 1; i < m; i++) {
            L->data[i*m+j] = LU.data[i*n+j];
        }
    }
    for (int i = 0; i < m; i++) { // Set diagonal elements of L to 1
        L->data[i*m+i] = 1;
    }

    for (int i = 0; i < rank; i++) {// Populate the upper triangular part of U with appropriate values from LU
        for (int j = 0; j < n; j++) {
            U->data[i*n+j] = LU.data[i*n+j];
        }
    }
    for (int i = 0; i < rank; i++) {// Set elements below the diagonal of U to 0
        for (int j = 0; j < i; j++) {
            U->data[i*n+j] = 0;
        }
    }
}

