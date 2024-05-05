#include "permutations.h"

// Swaps two rows of a matrix
void swapRows(Matrix* A, int row1, int row2) {
    int temp;
    for (int i = 0; i < A->cols; i++) {
        temp = A->data[row1*A->cols+i];
        A->data[row1*A->cols+i] = A->data[row2*A->cols+i];
        A->data[row2*A->cols+i] = temp;
    }
}

// Swaps two columns of a matrix
void swapCols(Matrix* A, int col1, int col2) {
    int temp;
    for (int i = 0; i < A->rows; i++) {
        temp = A->data[i*A->cols+col1];
        A->data[i*A->cols+col1] = A->data[i*A->cols+col2];
        A->data[i*A->cols+col2] = temp;
    }
}

// Permutes the rows of a matrix according to a given permutation array
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

// Permutes the columns of a matrix according to a given permutation array
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

// Performs a row transposition in a matrix, swapping two rows and updating a permutation array
void rowTransposition(Matrix* A,int src,int tgt,int* R){
    swapRows(A, src, tgt);
    int tmp = R[src];
    R[src] = R[tgt];
    R[tgt] = tmp;
}

// Performs a column transposition in a matrix, swapping two columns and updating a permutation array
void colTransposition(Matrix* A,int src,int tgt,int* C){
    swapCols(A, src, tgt);
    int tmp = C[src];
    C[src] = C[tgt];
    C[tgt] = tmp;
}

// Performs a row rotation in a matrix, shifting the rows and updating a permutation array ('row' to the last position)
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

// Performs a column rotation in a matrix, shifting the columns and updating a permutation array ('start_col' to 'end_col')
void colRotation(Matrix* A, int start_col, int end_col, int* R) {
    int cols = A->cols;
    int rows = A->rows;
    int* tempCol = (int*)malloc(rows * sizeof(int));
    // Copy the column to be rotated to a temporary array
    for (int i = 0; i < rows; i++) {
        tempCol[i] = A->data[i * cols + start_col];
    }
    int temp = R[start_col];
    //shift the columns
    for (int i = start_col; i < end_col-1; i++) {
        for (int row = 0; row < rows; row++)
            A->data[row * cols + i] = A->data[row * cols + (i + 1)];
        R[i] = R[i + 1];
    }
    // Copy the temporary array to 'end_col'
    for (int row = 0; row < rows; row++)
        A->data[row * cols + (end_col-1)] = tempCol[row];
    R[end_col-1] = temp;
    free(tempCol);
}
