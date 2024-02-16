#include "matrix.h"

Matrix createMatrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (int*)malloc(rows * cols * sizeof(int));
    return mat;
}

Matrix randomMatrix(int rows, int cols, int prime) {
    Matrix mat = createMatrix(rows, cols);
    srand(time(NULL));
    for (int i = 0; i < rows*cols; i++) {
            mat.data[i] = rand() % prime;
    }
    return mat;
}

void printMatrix(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            printf("%d\t", mat.data[i*mat.cols+j]);
        }
        printf("\n");
    }
}

void copyMatrix(Matrix src, Matrix* dest) {
    for (int i = 0; i < src.rows*src.cols; i++) {
            dest->data[i] = src.data[i];
    }
}

Matrix zerosMatrix(int rows, int cols) {
    Matrix mat = createMatrix(rows, cols);
    for (int i = 0; i < rows*cols; i++) {
            mat.data[i] = 0;
    }
    return mat;
}

Matrix multiplyMatrices(Matrix A, Matrix B, int prime) {
    if (A.cols != B.rows) {
        printf("Cannot multiply matrices: incompatible dimensions\n");
        exit(1);
    }

    Matrix result = createMatrix(A.rows, B.cols);

    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.cols; j++) {
            long long sum = 0;
            for (int k = 0; k < A.cols; k++) {
                sum += (long long) A.data[i*A.cols+k] * B.data[k*B.cols+j]; // 
            }
            result.data[i*B.cols+j] = sum % prime;
        }
    }

    return result;
}

bool compareMatrices(Matrix A, Matrix B) {
    if (A.rows != B.rows || A.cols != B.cols) {
        return false; // Matrices have different dimensions
    }

    for (int i = 0; i < A.rows*A.cols; i++) {
            if (A.data[i] != B.data[i]) {
                return false; // Matrices have different elements
            }
    }
    return true; // Matrices are equal
}
