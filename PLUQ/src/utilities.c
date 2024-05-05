#include "utilities.h"

// Returns the minimum of two integers
inline int min(int a, int b) { a < b ? a : b;}

// Allocates memory for an array of integers containing a range from 0 to (n - 1)
int* createRange(int n) {
    int* range = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        range[i] = i;
    }
    return range;
}

// Prints the elements of an integer array
void printArray(int* array, int size){
    for(int i = 0; i < size; i++){
        printf("%d\t", array[i]);
    }
    printf("\n");
}

// Checks if the given matrix is a lower triangular matrix
bool checkTriL(Matrix L) {
    int m = L.rows;
    int n = L.cols;
    if (m != n) {// If the matrix is not square
        return false;
    }
    for (int i = 0; i < m; i++) {
        if (L.data[i*n+i] != 1) {// Check if diagonal elements are 1
            return false;
        }
        for (int j = i + 1; j < m; j++) {// Check if elements above diagonal are 0
            if (L.data[i*n+j] != 0) {
                return false;
            }
        }
    }
    return true;
}

// Checks if the given matrix is an upper triangular matrix
bool checkTriU(Matrix U, int rank) {
    int m = U.rows;
    int n = U.cols;

    if (rank > m || rank > n) {// If rank exceeds matrix dimensions
        return false;
    }

    for (int i = 0; i < rank; i++) {
        if (U.data[i*n+i] == 0) {// Check if diagonal elements are nonzero
            return false;
        }
        for (int j = 0; j < min(i, n); j++) {
            if (U.data[i*n+j] != 0) {// Check if elements below diagonal are 0
                return false;
            }
        }
    }

    return true;
}

// Perform multiple tests on the PLUQ decomposition with different matrix sizes and permutations.
void checkManyPLUQ(int p, int max_iter, int algo){
    for (long test_case = 0; test_case < 3; test_case++)
    {
        int i = 0;
        bool correct = true;
        while (correct && i < max_iter){
            Matrix A;
            if (test_case == 0)
                A = randomMatrix(4, 4, p); // Generate random 4x4 matrix
            else if (test_case == 1)
                A = randomMatrix(3, 6, p); // Generate random 3x6 matrix
            else if (test_case == 2)
                A = randomMatrix(6, 3, p); // Generate random 6x3 matrix

            int* P = NULL;
            int* Q = NULL;
            Matrix LU;
            Matrix L;
            Matrix U;
            int rank = 0;
    
            PLUQ(A, &P, &LU, &Q, &rank, p, algo); // Perform PLUQ decomposition
            expand_PLUQ(LU, rank, &L, &U); // Expand LU matrix into L and U matrices
            correct = checkTriL(L) && checkTriU(U,rank); // Check if L is lower triangular and U is upper triangular
            permuteMatrixRows(&L, P); // Permute rows of L and columns of U according to permutation matrices P and Q
            permuteMatrixCols(&U, Q);
            Matrix L_mult_U = multiplyMatrices(L, U, p); // Multiply L and U matrices and compare the result with the original matrix A
            correct = correct && compareMatrices(L_mult_U, A); // check if A = LU
            i += 1;

            if (!correct)
            {
                printf("test_case %ld: wrong\n", test_case);
                return;
            }
        }
    }
    printf("---all tests passed---\n");
}

// Perform a single test on the PLUQ decomposition for a given matrix size and permutation.
void checkOnePLUQ(int p, int m, int n, bool print, int algo){
    bool correct = true;
    Matrix A = randomMatrix(m, n, p); // Generate random matrix A of size m x n
    int* P = NULL;
    int* Q = NULL;
    Matrix LU;
    Matrix L;
    Matrix U;
    int rank = 0;
    PLUQ(A, &P, &LU, &Q, &rank, p, algo); // Perform PLUQ decomposition
    expand_PLUQ(LU, rank, &L, &U); // Expand PLUQ into L and U matrices
        
    if(print){// Print matrices and permutations if print flag is true
        printf("A:\n");
        printMatrix(A);
        printf("P:\n");
        printArray(P, m);
        printf("L:\n");
        printMatrix(L);
        printf("U:\n");
        printMatrix(U);
        printf("Q:\n");
        printArray(Q, n);
        printf("\n");
    }

    correct = checkTriL(L) && checkTriU(U,rank); // Check if L is lower triangular and U is upper triangular
    // Permute rows of L and columns of U according to permutation matrices P and Q
    permuteMatrixRows(&L, P);
    permuteMatrixCols(&U, Q);
    // Multiply L and U matrices and compare the result with the original matrix A
    Matrix L_mult_U = multiplyMatrices(L, U, p);
    correct = correct && compareMatrices(L_mult_U, A);

    if (correct)
        printf("Matrix (%d, %d): ok\n", m, n);
    else{
        printf("Matrix (%d, %d): wrong\n", m, n);
    }

}
