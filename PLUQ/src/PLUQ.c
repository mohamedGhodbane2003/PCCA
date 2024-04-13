#include "PLUQ.h"

void pluq_crout(Matrix* A, int** P, int** Q, int* rank, int p) {
    int m = A->rows;
    int n = A->cols;

    *P = createRange(m);
    *Q = createRange(n);

    int matrixRank = 0;
    int nullity = 0;

    while(matrixRank + nullity < m) {
        for(int i = matrixRank; i < n; i++) {
            int tmp = 0;
            for(int j = 0; j < matrixRank; j++){
                tmp = add(tmp, mult(A->data[matrixRank * n + j], A->data[j * n + i], p), p);
            }
            A->data[matrixRank * n + i] = sub(A->data[matrixRank * n + i], tmp, p);
        }

        int pivot = matrixRank;
        while(pivot < n && A->data[matrixRank, pivot] == 0)
            pivot++;
        
        if(pivot == n) {
            rowRotation(A, matrixRank, *P);
            nullity++;
        }else{
            int invpivot = inverse(A->data[matrixRank, pivot], p);
            for(int i = matrixRank + 1; i < m - nullity; i++){
                int tmp = 0;
                 for(int j = 0; j < matrixRank; j++){
                    tmp = add(tmp, mult(A->data[i * n + j], A->data[j * n + pivot], p), p);
                 }
                 A->data[i * n + pivot] = sub(A->data[i * n + pivot], tmp, p);
                 A->data[i * n + pivot] = mult(A->data[i * n + pivot], invpivot, p);
            }
            colRotation(A, matrixRank, pivot+1, *Q);
            matrixRank++;
        }

    }
    *rank = matrixRank;
}

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
            __m256d vu = _mm256_set1_pd(1.0 / p);
            __m256d vp = _mm256_set1_pd(p);
            __m128i vp_128 = _mm_set1_epi32(p);
            for (int k = matrixRank + 1; k < m; k++) {
                A->data[k * n + matrixRank] = mult(A->data[k * n + matrixRank], inv, p);
                rows_elimination_avx2(A->data, n, matrixRank,A->data[k * n + matrixRank], p ,vp ,vu, vp_128 ,k);
            }
            matrixRank++;
        }
    }
    *rank = matrixRank;
}

void pluq_inplace(Matrix* A, int** P, int** Q, int* rank, int p) {
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
            for (int k = matrixRank + 1; k < m; k++) {
                A->data[k * n + matrixRank] = mult(A->data[k * n + matrixRank], inv, p);
                for (int j = matrixRank + 1; j < n; j++)
                    A->data[k * n + j] = sub(A->data[k * n + j], mult(A->data[k * n + matrixRank], A->data[matrixRank * n + j], p), p);
            }
            matrixRank++;
        }
    }
    *rank = matrixRank;
}

void PLUQ(Matrix A, int** P, Matrix* LU, int** Q, int* rank, int p) {
    *LU = createMatrix(A.rows, A.cols);
    copyMatrix(A, LU);
    //pluq_inplace(LU, P, Q, rank, p);
    //pluq_inplace_avx2(LU, P, Q, rank, p);
    pluq_crout(LU, P, Q, rank, p);
}

void expand_PLUQ(Matrix LU, int rank, Matrix* L, Matrix* U) {
    int m = LU.rows;
    int n = LU.cols;

    *L = zerosMatrix(m, m); 
    *U = zerosMatrix(m, n); 

    for (int j = 0; j < rank; j++) {
        for (int i = j + 1; i < m; i++) {
            L->data[i*m+j] = LU.data[i*n+j];
        }
    }
    for (int i = 0; i < m; i++) {
        L->data[i*m+i] = 1;
    }

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

