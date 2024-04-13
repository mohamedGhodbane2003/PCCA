#ifndef PERMUTATIONS_H
#define PERMUTATIONS_H

#include "matrix.h"

void swapRows(Matrix* A, int row1, int row2);
void swapCols(Matrix* A, int col1, int col2);
void permuteMatrixRows(Matrix* A, int* perm);
void permuteMatrixCols(Matrix* A, int* perm);
void rowTransposition(Matrix* A, int src, int tgt, int* R);
void colTransposition(Matrix* A, int src, int tgt, int* C);
void rowRotation(Matrix* A, int row, int* R);
void colRotation(Matrix* A, int start_col, int end_col, int* R);
#endif
