#include "permutations.h"



/*
 * Function: swapRows
 * -------------------
 * Swaps two rows of a matrix.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - row1: Index of the first row to swap.
 *   - row2: Index of the second row to swap.
 * 
 * Details:
 *   - Iterates through each element of the rows and swaps their values.
 */

void swapRows(Matrix* A, int row1, int row2) {
    int temp;
    for (int i = 0; i < A->cols; i++) {
        temp = A->data[row1][i];
        A->data[row1][i] = A->data[row2][i];
        A->data[row2][i] = temp;
    }
}

/*
 * Function: swapCols
 * -------------------
 * Swaps two columns of a matrix.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - col1: Index of the first column to swap.
 *   - col2: Index of the second column to swap.
 * 
 * Details:
 *   - Iterates through each element of the columns and swaps their values.
 */

void swapCols(Matrix* A, int col1, int col2) {
    int temp;
    for (int i = 0; i < A->rows; i++) {
        temp = A->data[i][col1];
        A->data[i][col1] = A->data[i][col2];
        A->data[i][col2] = temp;
    }
}

/*
 * Function: permuteMatrixRows
 * ----------------------------
 * Permutes the rows of a matrix based on given permutations.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - perm: Permutation array specifying the new order of rows.
 * 
 * Details:
 *   - Reorders the rows of the matrix according to the given permutation array.
 *   - Swaps rows until each row is in its correct position as per the permutation.
 */

void permuteMatrixRows(Matrix* A, int* perm) {
    for (int i = 0; i < A->rows; i++) {
        // If the current row is not already in its correct position
        while (perm[i] != i) {
            // Swap the current row with the row it should be in
            swapRows(A, i, perm[i]);
            // Update the permutation array to reflect the change
            int temp = perm[i];
            perm[i] = perm[temp];
            perm[temp] = temp;
        }
    }
}

/*
 * Function: permuteMatrixCols
 * ----------------------------
 * Permutes the columns of a matrix based on given permutations.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - perm: Permutation array specifying the new order of columns.
 * 
 * Details:
 *   - Reorders the columns of the matrix according to the given permutation array.
 *   - Swaps columns until each column is in its correct position as per the permutation.
 */

void permuteMatrixCols(Matrix* A, int* perm) {
    for (int i = 0; i < A->cols; i++) {
        // If the current column is not already in its correct position
        while (perm[i] != i) {
            // Swap the current column with the column it should be in
            swapCols(A, i, perm[i]);
            // Update the permutation array to reflect the change
            int temp = perm[i];
            perm[i] = perm[temp];
            perm[temp] = temp;
        }
    }
}

/*
 * Function: rowTransposition
 * ---------------------------
 * Performs a row transposition operation on a matrix.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - src: Index of the source row to swap.
 *   - tgt: Index of the target row to swap with the source row.
 *   - R: Permutation array for row transposition.
 * 
 * Details:
 *   - Swaps the specified source and target rows in the matrix.
 *   - Updates the permutation array accordingly.
 */

void rowTransposition(Matrix* A,int src,int tgt,int* R){
    swapRows(A, src, tgt);
    int tmp = R[src];
    R[src] = R[tgt];
    R[tgt] = tmp;
}

/*
 * Function: colTransposition
 * ---------------------------
 * Performs a column transposition operation on a matrix.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - src: Index of the source column to swap.
 *   - tgt: Index of the target column to swap with the source column.
 *   - C: Permutation array for column transposition.
 * 
 * Details:
 *   - Swaps the specified source and target columns in the matrix.
 *   - Updates the permutation array accordingly.
 */

void colTransposition(Matrix* A,int src,int tgt,int* C){
    swapCols(A, src, tgt);
    int tmp = C[src];
    C[src] = C[tgt];
    C[tgt] = tmp;
}

/*
 * Function: rowRotation
 * ---------------------
 * Rotates a row of a matrix by moving it to the last row position and shifting all rows below it up.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - row: Index of the row to rotate.
 *   - R: Permutation array for row rotation.
 * 
 * Details:
 *   - Allocates memory for a temporary array to hold the row being rotated.
 *   - Copies the row to be rotated to the temporary array.
 *   - Shifts all rows below the specified row one position up, wrapping the last row to the position of the rotated row.
 *   - Updates the permutation array accordingly to reflect the row rotation.
 */

void rowRotation(Matrix* A, int row, int* R) {

    int cols = A->cols;
    int rows = A->rows;

    int* tempRow = (int*)malloc(cols * sizeof(int));
    if (tempRow == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    // Copy the row to be rotated to a temporary array
    for (int i = 0; i < cols; i++) {
        tempRow[i] = A->data[row][i];
    }
    int temp = R[row];
    // Shift the rows
    for (int i = row; i < rows - 1; i++) {
        for(int col = 0; col < cols; col++)
            A->data[i][col] = A->data[i+1][col];
        R[i] = R[i + 1];
    }

    // Copy the original row to the last row
    for(int col = 0; col < cols; col++)
        A->data[rows - 1][col] = tempRow[col];
    R[rows - 1] = temp;

    free(tempRow);
}

/*
 * Function: colRotation
 * ---------------------
 * Rotates a column of a matrix by moving it to the last column position and shifting all columns to its right.
 * 
 * Parameters:
 *   - A: Pointer to the matrix.
 *   - col: Index of the column to rotate.
 *   - R: Permutation array for column rotation.
 * 
 * Details:
 *   - Allocates memory for a temporary array to hold the column being rotated.
 *   - Copies the column to be rotated to the temporary array.
 *   - Shifts all columns to the right of the specified column one position left, wrapping the last column to the position of the rotated column.
 *   - Updates the permutation array accordingly to reflect the column rotation.
 */

void colRotation(Matrix* A, int col, int* R) {
    int cols = A->cols;
    int rows = A->rows;
    int* tempCol = (int*)malloc(rows * sizeof(int));
    if (tempCol == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }

    // Copy the column to be rotated to a temporary array
    for (int i = 0; i < rows; i++) {
        tempCol[i] = A->data[i][col];
    }
    int temp = R[col];
    // Shift the rows
    for (int i = col; i < cols - 1; i++) {
        for(int row = 0; row < rows; row++)
            A->data[row][i] = A->data[row][i+1];
        R[i] = R[i + 1];
    }

    // Copy the original row to the last row
    for(int row = 0; row < rows; row++)
        A->data[row][cols - 1] = tempCol[row];
    R[cols - 1] = temp;

    free(tempCol);
}