"""
In the Gaussian Elimination algorithm, 
the forward elimination phase aims to transform the matrix into an upper triangular form, 
performing row operations to transform the matrix into echelon form, 
This includes finding a nonzero pivot element in the current column, 
swapping rows to place the pivot element in the correct position, 
normalizing the pivot row, and eliminating entries below the pivot.
"""
def gaussian_elimination(matrix):
    """
    Input: The function takes a matrix as its input parameter
    Output: It returns the echelon form of the input matrix
    """

    rows, cols = matrix.nrows(), matrix.ncols()

    # Create a new matrix over the same base ring as the input matrix
    echelon_form = Matrix(matrix.base_ring(), rows, cols)

    # Copy the entries to the new matrix
    for i in range(rows):
        for j in range(cols):
            echelon_form[i, j] = matrix[i, j]

    # Forward Elimination (using min(rows, cols - 1) ensures that it won't attempt 
    # to access rows beyond the actual number of rows in the matrix)

    for col in range(min(rows, cols - 1)):
        pivot_row = -1

        # Find the pivot row with maximum absolute value in the column
        max_abs_value = 0
        for i in range(col, rows):
            if abs(echelon_form[i, col]) > max_abs_value:
                max_abs_value = abs(echelon_form[i, col])
                pivot_row = i

        # If no nonzero pivot is found, move to the next column
        if pivot_row == -1:
            continue

        # Swap rows to move the pivot to the current row
        echelon_form.swap_rows(col, pivot_row)

        pivot_row = col

        # Eliminate entries below the pivot
        for i in range(pivot_row + 1, rows):
            factor = echelon_form[i, col] / echelon_form[pivot_row, col]
            echelon_form[i, :] -= factor * echelon_form[pivot_row, :]

    return echelon_form



# Example Usage
# Create a matrix over the Rational Field
A = Matrix(QQ, [
    [ 1/2,   -1, -1/2,  1/2],
    [  -2,    0,    0,    0],
    [   2,    0,    1,  1/2],
    [  -1, -1/2,   -1,   -1]
    
])



# Apply Gaussian Elimination
echelon_form_A = gaussian_elimination(A)
P, L, U = A.LU()

# Display the result
print("Original Matrix A:")
print(A)
print("\nEchelon Form of A:")
print(echelon_form_A)
print("\nU matrix in the LU(A):")
print(U)


# example over finite field Z/pZ:
# A = matrix.random(GF(3), 4, 4)
