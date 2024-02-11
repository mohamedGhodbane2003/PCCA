#include "PLUQ.h"

/*
 * Function: PLUQ
 * --------------
 * Performs PLUQ factorization of a matrix A over a finite field.
 * 
 * Parameters:
 *   - A: Pointer to the input matrix to be factorized.
 *   - P: Pointer to the permutation matrix P.
 *   - LU: Pointer to the matrix LU (L-U decomposition of A).
 *   - Q: Pointer to the permutation matrix Q.
 *   - rank: Pointer to an integer to store the rank of the matrix A.
 *   - p: Prime number used for modulo operations.
 * 
 * Details:
 *   PLUQ factorization decomposes the input matrix A into three matrices:
 *   - P: A permutation matrix that represents the row permutations performed during factorization.
 *   - LU: A matrix containing both the lower triangular matrix (L) and upper triangular matrix (U).
 *   - Q: A permutation matrix that represents the column permutations performed during factorization.
 * 
 *   The decomposition is such that A = PLUQ.
 * 
 *   The algorithm iteratively transforms A into an upper triangular matrix (U) using Gaussian elimination.
 *   At each step, it also keeps track of row and column permutations (P and Q) and computes the lower triangular matrix (L).
 * 
 *   The factorization is performed over a finite field with modulus p, ensuring that all intermediate computations
 *   are within the finite field.
 * 
 *   The rank of the matrix A is determined during the factorization process and stored in the 'rank' variable.
 * 
 *   Steps:
 *   - Initialize LU with a copy of A.
 *   - Initialize permutations P and Q with identity permutations.
 *   - Perform Gaussian elimination to transform LU into an upper triangular matrix (U).
 *   - During Gaussian elimination:
 *     - Find the pivot element in the current column by searching for the first non-zero element below the current row.
 *     - If no pivot is found, perform row rotation and increase the nullity count.
 *     - If a pivot is found:
 *       - Perform column transposition to move the pivot element to the diagonal position.
 *       - Use Gaussian elimination to eliminate all elements below the pivot.
 *       - Adjust the matrix rank and continue to the next column.
 *   - After completing Gaussian elimination, the resulting LU matrix represents both L and U.
 *   - The rank of A is determined by counting the number of non-zero rows in U.
 * 
 */


void PLUQ(Matrix* A, int** P, Matrix** LU, int** Q, int* rank, int p){
    int m = A->rows;
    int n = A->cols;
    *LU = createMatrix(m, n);
    // Initialize LU with copy of A
    copyMatrix(A, *LU);
    // Initialize permutations with identity permutation
    *Q = createRange(n);
    *P = createRange(m);

    // Currently discovered rank and dimension of left nullspace
    int matrixRank = 0;
    int nullity = 0;

    while (matrixRank + nullity < m){
        // Find column with pivot element on row `matrixRank`, if there is some
        int pivot = matrixRank; // pivot is at column index >= matrixRank; take first one
        while (pivot < n && (*LU)->data[matrixRank][pivot] == 0)
            pivot += 1;
        if (pivot == n){
            rowRotation(*LU,matrixRank,*P);
            nullity++;
         }else{
            colTransposition(*LU,matrixRank,pivot,*Q);
            for (int k = matrixRank+1; k < m ; k++){
                (*LU)->data[k][matrixRank] = mult((*LU)->data[k][matrixRank], inverse((*LU)->data[matrixRank][matrixRank], p), p);
                for(int j = matrixRank+1; j < n; j++)
                    (*LU)->data[k][j] = sub((*LU)->data[k][j], mult((*LU)->data[k][matrixRank], (*LU)->data[matrixRank][j], p), p);
            }
            matrixRank++;
         }
    }
    *rank = matrixRank;
}



/*
 * Function: expand_PLUQ
 * ----------------------
 * Expands the LU decomposition into separate matrices L and U.
 * 
 * Parameters:
 *   - LU: Pointer to the LU decomposition matrix (combined L and U).
 *   - rank: Rank of the original matrix A.
 *   - L: Pointer to store the resulting lower triangular matrix L.
 *   - U: Pointer to store the resulting upper triangular matrix U.
 * 
 * Details:
 *   Given the LU decomposition matrix (LU), which represents both L and U, this function
 *   separates it into two matrices: L (lower triangular) and U (upper triangular).
 * 
 *   The dimensions of L are m x m, where m is the number of rows in LU, and U has the same
 *   dimensions as LU.
 * 
 *   The expansion is based on the rank of the original matrix A. Only the elements relevant
 *   to the rank are retained in both L and U.
 * 
 *   After expansion, L will be a proper lower triangular matrix, and U will be an upper
 *   triangular matrix with zeros below the diagonal.
 * 
 */


void expand_PLUQ(Matrix* LU, int rank, Matrix** L, Matrix** U) {
    int m = LU->rows;
    int n = LU->cols;

    *L = zerosMatrix(m, m); // L is square
    *U = zerosMatrix(m, n); // U has the same dimensions as LU

    // Retrieve L
    for (int j = 0; j < rank; j++) {
        for (int i = j + 1; i < m; i++) {
            (*L)->data[i][j] = LU->data[i][j];
        }
    }
    for (int i = 0; i < m; i++) {
        (*L)->data[i][i] = 1;
    }

    // Retrieve U
    for (int i = 0; i < rank; i++) {
        for (int j = 0; j < n; j++) {
            (*U)->data[i][j] = LU->data[i][j];
        }
    }
    for (int i = 0; i < rank; i++) {
        for (int j = 0; j < i; j++) {
            (*U)->data[i][j] = 0;
        }
    }
}