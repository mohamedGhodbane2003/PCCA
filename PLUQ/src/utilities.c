#include "utilities.h"

inline int min(int a, int b) { a < b ? a : b;}

int* createRange(int n) {
    int* range = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        range[i] = i;
    }
    return range;
}

void printArray(int* array, int size){
    for(int i = 0; i < size; i++){
        printf("%d\t", array[i]);
    }
    printf("\n");
}

bool checkTriL(Matrix L) {
    int m = L.rows;
    int n = L.cols;
    if (m != n) {
        return false;
    }
    for (int i = 0; i < m; i++) {
        if (L.data[i*n+i] != 1) {
            return false;
        }
        for (int j = i + 1; j < m; j++) {
            if (L.data[i*n+j] != 0) {
                return false;
            }
        }
    }
    return true;
}

bool checkTriU(Matrix U, int rank) {
    int m = U.rows;
    int n = U.cols;

    if (rank > m || rank > n) {
        return false;
    }

    for (int i = 0; i < rank; i++) {
        if (U.data[i*n+i] == 0) {
            return false;
        }
        for (int j = 0; j < min(i, n); j++) {
            if (U.data[i*n+j] != 0) {
                return false;
            }
        }
    }

    return true;
}

void checkManyPLUQ(int p, int max_iter){
    for (long test_case = 0; test_case < 3; test_case++)
    {
        int i = 0;
        bool correct = true;
        while (correct && i < max_iter){
            Matrix A;
            if (test_case == 0)
                A = randomMatrix(4, 4, p);
            else if (test_case == 1)
                A = randomMatrix(3, 6, p);
            else if (test_case == 2)
                A = randomMatrix(6, 3, p);

            int* P = NULL;
            int* Q = NULL;
            Matrix LU;
            Matrix L;
            Matrix U;
            int rank = 0;

            PLUQ(A, &P, &LU, &Q, &rank, p);
            expand_PLUQ(LU, rank, &L, &U);
            correct = checkTriL(L) && checkTriU(U,rank);
            permuteMatrixRows(&L, P);
            permuteMatrixCols(&U, Q);
            Matrix L_mult_U = multiplyMatrices(L, U, p);
            correct = correct && compareMatrices(L_mult_U, A);
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

void checkOnePLUQ(int p, int m, int n, bool print){
    
    bool correct = true;
    Matrix A = randomMatrix(m, n, p);
    int* P = NULL;
    int* Q = NULL;
    Matrix LU;
    Matrix L;
    Matrix U;
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
    permuteMatrixRows(&L, P);
    permuteMatrixCols(&U, Q);
    Matrix L_mult_U = multiplyMatrices(L, U, p);
    correct = correct && compareMatrices(L_mult_U, A);

    if (correct)
        printf("Matrix (%d, %d): ok\n", m, n);
    else{
        printf("Matrix (%d, %d): wrong\n", m, n);
    }

}
