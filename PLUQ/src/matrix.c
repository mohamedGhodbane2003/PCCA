#include "matrix.h"

/*
 * Function: createMatrix
 * ----------------------
 * Creates a matrix with the specified number of rows and columns.
 * 
 * Parameters:
 *   - rows: The number of rows in the matrix.
 *   - cols: The number of columns in the matrix.
 * 
 * Returns:
 *   - A pointer to the created matrix.
 * 
 * Memory Allocation:
 *   - Allocates memory for the Matrix struct.
 *   - Allocates memory for the data array, representing the matrix elements.
 *   - Allocates memory for each row in the data array.
 * 
 * Error Handling:
 *   - Prints an error message and exits the program if memory allocation fails.
 */

Matrix* createMatrix(int rows, int cols) {
    Matrix* mat = (Matrix*)malloc(sizeof(Matrix));
    if (mat == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    mat->rows = rows;
    mat->cols = cols;
    mat->data = (int**)malloc(rows * sizeof(int*));
    if (mat->data == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < rows; i++) {
        mat->data[i] = (int*)malloc(cols * sizeof(int));
        if (mat->data[i] == NULL) {
            printf("Memory allocation failed\n");
            exit(1);
        }
    }
    return mat;
}

/*
 * Function: randomMatrix
 * ----------------------
 * Creates a matrix with the specified number of rows and columns, filled with random numbers 
 * in the range [0, prime - 1], where prime is a given parameter.
 * 
 * Parameters:
 *   - rows: The number of rows in the matrix.
 *   - cols: The number of columns in the matrix.
 *   - prime: The upper limit for the random numbers (exclusive).
 * 
 * Returns:
 *   - A pointer to the created matrix.
 * 
 * Memory Allocation:
 *   - Calls the createMatrix function to allocate memory for the matrix.
 * 
 * Random Number Generation:
 *   - Seeds the random number generator using the current time.
 *   - Fills each element of the matrix with a random number in the range [0, prime - 1].
 * 
 */

Matrix* randomMatrix(int rows, int cols, int prime) {
    Matrix* mat = createMatrix(rows, cols);
    srand(time(NULL));
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] = rand() % prime;
        }
    }
    return mat;
}


/*
 * Function: destroyMatrix
 * -----------------------
 * Deallocates the memory used by a matrix, including the matrix itself and its data array.
 * 
 * Parameters:
 *   - mat: A pointer to the matrix to be destroyed.
 * 
 * Memory Deallocation:
 *   - Frees the memory allocated for each row in the data array.
 *   - Frees the memory allocated for the data array itself.
 *   - Frees the memory allocated for the matrix structure.
 * 
 */

void destroyMatrix(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        free(mat->data[i]);
    }
    free(mat->data);
    free(mat);
}

/*
 * Function: printMatrix
 * ----------------------
 * Prints the elements of a matrix to the standard output.
 * 
 * Parameters:
 *   - mat: A pointer to the matrix to be printed.
 * 
 * Printing Format:
 *   - The elements of the matrix are printed row by row, separated by tabs (\t).
 *   - Each row is printed on a separate line.
 * 
 * Example:
 *   If mat is a 2x3 matrix with elements:
 *   1 2 3
 *   4 5 6
 *   
 *   The output of printMatrix(mat) will be:
 *   1  2  3
 *   4  5  6
 */

void printMatrix(Matrix* mat) {
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            printf("%d\t", mat->data[i][j]);
        }
        printf("\n");
    }
}

/*
 * Function: copyMatrix
 * ---------------------
 * Copies the contents of one matrix to another matrix.
 * 
 * Parameters:
 *   - src: A pointer to the source matrix from which data is copied.
 *   - dest: A pointer to the destination matrix where data is copied to.
 * 
 * Notes:
 *   - The dimensions (rows and columns) of the source and destination matrices must match.
 *   - The function performs a deep copy, copying each individual element from src to dest.
 * 
 * Example:
 *   If src is a 2x3 matrix with elements:
 *   1 2 3
 *   4 5 6
 *   
 *   After calling copyMatrix(src, dest), where dest is an initially empty matrix,
 *   dest will also become a 2x3 matrix with elements:
 *   1 2 3
 *   4 5 6
 */

void copyMatrix(Matrix* src, Matrix* dest) {
    for (int i = 0; i < src->rows; i++) {
        for (int j = 0; j < src->cols; j++) {
            dest->data[i][j] = src->data[i][j];
        }
    }
}

/*
 * Function: zerosMatrix
 * ----------------------
 * Creates a matrix with the specified number of rows and columns, filled with zeros.
 * 
 * Parameters:
 *   - rows: The number of rows in the matrix.
 *   - cols: The number of columns in the matrix.
 * 
 * Returns:
 *   A pointer to the newly created matrix filled with zeros.
 * 
 * 
 * Example:
 *   If rows = 2 and cols = 3, the function will create a 2x3 matrix with elements:
 *   0 0 0
 *   0 0 0
 */

Matrix* zerosMatrix(int rows, int cols) {
    Matrix* mat = (Matrix*)malloc(sizeof(Matrix));
    if (mat == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    mat->rows = rows;
    mat->cols = cols;
    mat->data = (int**)malloc(rows * sizeof(int*));
    if (mat->data == NULL) {
        printf("Memory allocation failed\n");
        exit(1);
    }
    for (int i = 0; i < rows; i++) {
        mat->data[i] = (int*)malloc(cols * sizeof(int));
        if (mat->data[i] == NULL) {
            printf("Memory allocation failed\n");
            exit(1);
        }
        // Initialize each element to zero
        for (int j = 0; j < cols; j++) {
            mat->data[i][j] = 0;
        }
    }
    return mat;
}

/*
 * Function: multiplyMatrices
 * ---------------------------
 * Multiplies two matrices and returns the result.
 * 
 * Parameters:
 *   - A: Pointer to the first matrix (left operand).
 *   - B: Pointer to the second matrix (right operand).
 *   - prime: Prime number used for modulo operation during multiplication.
 * 
 * Returns:
 *   A pointer to the resulting matrix of the multiplication operation.
 */

Matrix* multiplyMatrices(Matrix* A, Matrix* B, int prime) {
    if (A->cols != B->rows) {
        printf("Cannot multiply matrices: incompatible dimensions\n");
        exit(1);
    }

    int m = A->rows;
    int n = A->cols;
    int p = B->cols;

    Matrix* result = createMatrix(m, p);

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            long long sum = 0;
            for (int k = 0; k < n; k++) {
                sum += (long long) A->data[i][k] * B->data[k][j]; // 
            }
            result->data[i][j] = sum % prime;
        }
    }

    return result;
}

/*
 * Function: compareMatrices
 * --------------------------
 * Compares two matrices for equality.
 * 
 * Parameters:
 *   - A: Pointer to the first matrix.
 *   - B: Pointer to the second matrix.
 * 
 * Returns:
 *   - true if the matrices are equal (same dimensions and same elements), otherwise false.
 * 
 */

bool compareMatrices(Matrix* A, Matrix* B) {
    if (A->rows != B->rows || A->cols != B->cols) {
        return false; // Matrices have different dimensions
    }

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->cols; j++) {
            if (A->data[i][j] != B->data[i][j]) {
                return false; // Matrices have different elements
            }
        }
    }

    return true; // Matrices are equal
}