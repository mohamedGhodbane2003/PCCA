#include "PLUQ.h"

void PLUQ(Matrix A, int** P, Matrix* LU, int** Q, int* rank, int p){
    int m = A.rows;
    int n = A.cols;
    // Initialize permutations with identity permutation
    *Q = createRange(n);
    *P = createRange(m);

    // Currently discovered rank and dimension of left nullspace
    int matrixRank = 0;
    int nullity = 0;

    while (matrixRank + nullity < m){
        int pivot = matrixRank; // pivot is at column index >= matrixRank; take first one
        while (pivot < n && A.data[matrixRank*n+pivot] == 0)
            pivot += 1;
        if (pivot == n){
            rowRotation(&A,matrixRank,*P);
            nullity++;
         }else{
            int inv = inverse(A.data[matrixRank*n+matrixRank], p);
            if(pivot != matrixRank)
                colTransposition(&A,matrixRank,pivot,*Q);
            for (int k = matrixRank+1; k < m ; k++){
                A.data[k*n+matrixRank] = mult(A.data[k*n+matrixRank], inv, p);
                for(int j = matrixRank+1; j < n; j++)
                    A.data[k*n+j] = sub(A.data[k*n+j], mult(A.data[k*n+matrixRank], A.data[matrixRank*n+j], p), p);
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

