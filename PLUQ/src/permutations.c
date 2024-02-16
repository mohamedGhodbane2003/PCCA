#include "permutations.h"

void swapRows(Matrix* A, int row1, int row2) {
    int temp;
    for (int i = 0; i < A->cols; i++) {
        temp = A->data[row1*A->cols+i];
        A->data[row1*A->cols+i] = A->data[row2*A->cols+i];
        A->data[row2*A->cols+i] = temp;
    }
}

void swapCols(Matrix* A, int col1, int col2) {
    int temp;
    for (int i = 0; i < A->rows; i++) {
        temp = A->data[i*A->cols+col1];
        A->data[i*A->cols+col1] = A->data[i*A->cols+col2];
        A->data[i*A->cols+col2] = temp;
    }
}

void permuteMatrixRows(Matrix* A, int* perm) {
    for (int i = 0; i < A->rows; i++) {
        while (perm[i] != i) {
            swapRows(A, i, perm[i]);
            int temp = perm[i];
            perm[i] = perm[temp];
            perm[temp] = temp;
        }
    }
}

void permuteMatrixCols(Matrix* A, int* perm) {
    for (int i = 0; i < A->cols; i++) {
        while (perm[i] != i) {
            swapCols(A, i, perm[i]);
            int temp = perm[i];
            perm[i] = perm[temp];
            perm[temp] = temp;
        }
    }
}

void rowTransposition(Matrix* A,int src,int tgt,int* R){
    swapRows(A, src, tgt);
    int tmp = R[src];
    R[src] = R[tgt];
    R[tgt] = tmp;
}

void colTransposition(Matrix* A,int src,int tgt,int* C){
    swapCols(A, src, tgt);
    int tmp = C[src];
    C[src] = C[tgt];
    C[tgt] = tmp;
}

void rowRotation(Matrix* A, int row, int* R) {

    int cols = A->cols;
    int rows = A->rows;

    int* tempRow = (int*)malloc(cols * sizeof(int));

    // Copy the row to be rotated to a temporary array
    for (int i = 0; i < cols; i++) {
        tempRow[i] = A->data[row*cols+i];
    }
    int temp = R[row];
    // Shift the rows
    for (int i = row; i < rows - 1; i++) {
        for(int col = 0; col < cols; col++)
            A->data[i*cols+col] = A->data[(i+1)*cols+col];
        R[i] = R[i + 1];
    }

    // Copy the original row to the last row
    for(int col = 0; col < cols; col++)
        A->data[(rows-1)*cols+col] = tempRow[col];
    R[rows-1] = temp;

    free(tempRow);
}

void colRotation(Matrix* A, int col, int* R) {
    int cols = A->cols;
    int rows = A->rows;
    int* tempCol = (int*)malloc(rows * sizeof(int));

    // Copy the column to be rotated to a temporary array
    for (int i = 0; i < rows; i++) {
        tempCol[i] = A->data[i*cols+col];
    }
    int temp = R[col];
    // Shift the rows
    for (int i = col; i < cols - 1; i++) {
        for(int row = 0; row < rows; row++)
            A->data[row*cols+i] = A->data[row*cols+(i+1)];
        R[i] = R[i+1];
    }

    // Copy the original row to the last row
    for(int row = 0; row < rows; row++)
        A->data[row*cols+(cols-1)] = tempCol[row];
    R[cols-1] = temp;

    free(tempCol);
}