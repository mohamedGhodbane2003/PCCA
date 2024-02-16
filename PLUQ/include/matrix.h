#ifndef MATRIX_H
#define MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>

typedef struct {
    int rows;
    int cols;
    int* data;
} Matrix;

Matrix createMatrix(int rows, int cols);
Matrix randomMatrix(int rows, int cols, int prime);
void printMatrix(Matrix mat);
void copyMatrix(Matrix src, Matrix* dest);
Matrix zerosMatrix(int rows, int cols);
Matrix multiplyMatrices(Matrix A, Matrix B, int prime);
bool compareMatrices(Matrix A, Matrix B);

#endif
