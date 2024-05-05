#include "matrix.h"

// Creates a matrix with the specified number of rows and columns, and allocates memory for its data
Matrix createMatrix(int rows, int cols) {
    Matrix mat;
    mat.rows = rows;
    mat.cols = cols;
    mat.data = (int*)malloc(rows * cols * sizeof(int));
    return mat;
}

// Creates a matrix with the specified number of rows and columns, and fills it with random values modulo 'prime'
Matrix randomMatrix(int rows, int cols, int prime) {
    Matrix mat = createMatrix(rows, cols);
    srand(time(NULL));
    for (int i = 0; i < rows*cols; i++) {
            mat.data[i] = rand() % prime;
    }
    return mat;
}

// Prints the contents of a matrix
void printMatrix(Matrix mat) {
    for (int i = 0; i < mat.rows; i++) {
        for (int j = 0; j < mat.cols; j++) {
            printf("%d\t", mat.data[i*mat.cols+j]);
        }
        printf("\n");
    }
}

// Copies the contents of one matrix to another
void copyMatrix(Matrix src, Matrix* dest) {
    for (int i = 0; i < src.rows*src.cols; i++) {
            dest->data[i] = src.data[i];
    }
}

// Creates a matrix with the specified number of rows and columns, filled with zeros.
Matrix zerosMatrix(int rows, int cols) {
    Matrix mat = createMatrix(rows, cols);
    for (int i = 0; i < rows*cols; i++) {
            mat.data[i] = 0;
    }
    return mat;
}

// Multiplies two matrices 'A' and 'B' modulo 'prime' (used for testing)
Matrix multiplyMatrices(Matrix A, Matrix B, int prime) {
    Matrix result = createMatrix(A.rows, B.cols);
    for (int i = 0; i < A.rows; i++) {
        for (int j = 0; j < B.cols; j++) {
            long long sum = 0;
            for (int k = 0; k < A.cols; k++) {
                sum += (long long) A.data[i*A.cols+k] * B.data[k*B.cols+j];
            }
            result.data[i*B.cols+j] = sum % prime;
        }
    }
    return result;
}

// Compares two matrices 'A' and 'B' and returns true if they are equal, false otherwise.
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
