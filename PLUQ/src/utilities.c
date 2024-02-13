#include "utilities.h"

/*
 * Function: createRange
 * ----------------------
 * Creates an array of integers from 0 to n-1.
 * 
 * Parameters:
 *   - n: Number of integers to generate.
 * 
 * Returns:
 *   A pointer to the array of integers.
 */

int* createRange(int n) {
    int* range = (int*)malloc(n * sizeof(int));
    if (range == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        range[i] = i;
    }
    return range;
}

/*
 * Function: printArray
 * ---------------------
 * Prints the elements of an integer array.
 * 
 * Parameters:
 *   - array: Pointer to the integer array.
 *   - size: Size of the array.
 */

void printArray(int* array, int size){
    for(int i = 0; i < size; i++){
        printf("%d\t", array[i]);
    }
    printf("\n");
}

/*
 * Function: checkTriL
 * --------------------
 * Checks if the given matrix is lower triangular.
 * 
 * Parameters:
 *   - L: Pointer to the matrix to be checked.
 * 
 * Returns:
 *   True if the matrix is lower triangular, otherwise false.
 */

bool checkTriL(Matrix* L) {
    int m = L->rows;
    int n = L->cols;
    if (m != n) {
        return false;
    }
    for (int i = 0; i < m; i++) {
        if (L->data[i][i] != 1) {
            return false;
        }
        for (int j = i + 1; j < m; j++) {
            if (L->data[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

/*
 * Function: min
 * --------------
 * Returns the minimum of two integers.
 */

int min(int a, int b){
    return a < b ? a : b;
}

/*
 * Function: checkTriU
 * --------------------
 * Checks if the given matrix is upper triangular.
 * 
 * Parameters:
 *   - U: Pointer to the matrix to be checked.
 *   - rank: Rank of the matrix.
 * 
 * Returns:
 *   True if the matrix is upper triangular, otherwise false.
 */

bool checkTriU(Matrix* U, int rank) {
    int m = U->rows;
    int n = U->cols;
    if (rank > m || rank > n) {
        return false;
    }
    for (int i = 0; i < rank; i++) {
        if (U->data[i][i] == 0) {
            return false;
        }
        for (int j = 0; j < min(i, n); j++) {
            if (U->data[i][j] != 0) {
                return false;
            }
        }
    }
    return true;
}

/*
 * Function: checkManyPLUQ
 * ------------------------
 * Checks the correctness of PLUQ decomposition for different types of matrices.
 * 
 * Parameters:
 *   - p: Prime number used in the finite field.
 *   - max_iter: Maximum number of iterations to perform the check.
 * 
 * Details:
 *   - Generates random matrices of different shapes and sizes.
 *   - Computes the PLUQ decomposition for each matrix and checks if the decomposition
 *     is correct by reconstructing the original matrix.
 *   - Prints the result for each matrix type (tall rectangular, wide rectangular, square).
 */

void checkManyPLUQ(int p, int max_iter){
    int i = 0;
    bool correct = true;
    while (correct && i < max_iter){
        Matrix* A = randomMatrix(6, 3, p);
        int* P = NULL;
        int* Q = NULL;
        Matrix* LU = NULL;
        Matrix* L = NULL;
        Matrix* U = NULL;
        int rank = 0;
        PLUQ(A, &P, &LU, &Q, &rank, p);
        expand_PLUQ(LU, rank, &L, &U);
        correct = checkTriL(L) && checkTriU(U,rank);
        permuteMatrixRows(L, P);
        permuteMatrixCols(U, Q);
        Matrix* L_mult_U = multiplyMatrices(L, U, p);
        correct = correct && compareMatrices(L_mult_U, A);
        i += 1;
    }

    if (correct)
        printf("tall rectangular: ok\n");
    else{
        printf("tall rectangular: wrong\n");
    }

    i = 0;
    correct = true;
    while (correct && i < max_iter){
        Matrix* A = randomMatrix(3, 6, p);
        int* P = NULL;
        int* Q = NULL;
        Matrix* LU = NULL;
        Matrix* L = NULL;
        Matrix* U = NULL;
        int rank = 0;
        PLUQ(A, &P, &LU, &Q, &rank, p);
        expand_PLUQ(LU, rank, &L, &U);
        correct = checkTriL(L) && checkTriU(U,rank);
        permuteMatrixRows(L, P);
        permuteMatrixCols(U, Q);
        Matrix* L_mult_U = multiplyMatrices(L, U, p);
        correct = correct && compareMatrices(L_mult_U, A);
        i += 1;
    }

    if (correct)
        printf("wide rectangular: ok\n");
    else{
        printf("wide rectangular: wrong\n");
    }

    i = 0;
    correct = true;
    while (correct && i < max_iter){
        Matrix* A = randomMatrix(4, 4, p);
        int* P = NULL;
        int* Q = NULL;
        Matrix* LU = NULL;
        Matrix* L = NULL;
        Matrix* U = NULL;
        int rank = 0;
        PLUQ(A, &P, &LU, &Q, &rank, p);
        expand_PLUQ(LU, rank, &L, &U);
        correct = checkTriL(L) && checkTriU(U,rank);
        permuteMatrixRows(L, P);
        permuteMatrixCols(U, Q);
        Matrix* L_mult_U = multiplyMatrices(L, U, p);
        correct = correct && compareMatrices(L_mult_U, A);
        i += 1;
    }

    if (correct)
        printf("square: ok\n");
    else{
        printf("square: wrong\n");
    }

}

/*
 * Function: checkOnePLUQ
 * ------------------------
 * Checks the correctness of PLUQ decomposition for one matrix.
 * 
 * Parameters:
 *   - p: Prime number used in the finite field.
 *   - m: number of rows
 *   - n: number of columns.
 *   - print: flag to print the matrices
 * 
 * Details:
 *   - Generates random matrix (m, n)
 *   - Computes the PLUQ decomposition and checks if the decomposition
 *     is correct by reconstructing the original matrix.
 *   - Print the result
 */

void checkOnePLUQ(int p, int m, int n, bool print){
    
    bool correct = true;
    Matrix* A = randomMatrix(m, n, p);
    int* P = NULL;
    int* Q = NULL;
    Matrix* LU = NULL;
    Matrix* L = NULL;
    Matrix* U = NULL;
    int rank = 0;
    PLUQ(A, &P, &LU, &Q, &rank, p);
    expand_PLUQ(LU, rank, &L, &U);

    if(print){
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

    correct = checkTriL(L) && checkTriU(U,rank);
    permuteMatrixRows(L, P);
    permuteMatrixCols(U, Q);
    Matrix* L_mult_U = multiplyMatrices(L, U, p);
    correct = correct && compareMatrices(L_mult_U, A);

    if (correct)
        printf("Matrix (%d, %d): ok\n", m, n);
    else{
        printf("Matrix (%d, %d): wrong\n");
    }

}