#include "PLUQ.h"

void PLUQ(Matrix A, int** P, Matrix* LU, int** Q, int* rank, int p){
    int m = A.rows;
    int n = A.cols;
    *LU = createMatrix(m, n);
    // Initialize LU with copy of A
    copyMatrix(A, LU);

    // Initialize permutations with identity permutation
    *Q = createRange(n);
    *P = createRange(m);

    // Currently discovered rank and dimension of left nullspace
    int matrixRank = 0;
    int nullity = 0;

    while (matrixRank + nullity < m){
        int pivot = matrixRank; // pivot is at column index >= matrixRank; take first one
        while (pivot < n && LU->data[matrixRank*n+pivot] == 0)
            pivot += 1;
        if (pivot == n){
            rowRotation(LU,matrixRank,*P);
            nullity++;
         }else{
            int inv = inverse(LU->data[matrixRank*n+matrixRank], p);
            if(pivot != matrixRank)
                colTransposition(LU,matrixRank,pivot,*Q);
            for (int k = matrixRank+1; k < m ; k++){
                LU->data[k*n+matrixRank] = mult(LU->data[k*n+matrixRank], inv, p);
                for(int j = matrixRank+1; j < n; j++)
                    LU->data[k*n+j] = sub(LU->data[k*n+j], mult(LU->data[k*n+matrixRank], LU->data[matrixRank*n+j], p), p);
            }
            matrixRank++;
         }
    }
    *rank = matrixRank;
}

void expand_PLUQ(Matrix LU, int rank, Matrix* L, Matrix* U) {
    int m = LU.rows;
    int n = LU.cols;

    *L = zerosMatrix(m, m); // L is square
    *U = zerosMatrix(m, n); // U has the same dimensions as LU

    // Retrieve L
    for (int j = 0; j < rank; j++) {
        for (int i = j + 1; i < m; i++) {
            L->data[i*m+j] = LU.data[i*n+j];
        }
    }
    for (int i = 0; i < m; i++) {
        L->data[i*m+i] = 1;
    }

    // Retrieve U
    for (int i = 0; i < rank; i++) {
        for (int j = 0; j < n; j++) {
            U->data[i*n+j] = LU.data[i*n+j];
        }
    }
    for (int i = 0; i < rank; i++) {
        for (int j = 0; j < i; j++) {
            U->data[i*n+j] = 0;
        }
    }
}

