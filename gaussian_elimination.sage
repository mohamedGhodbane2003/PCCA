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

        # Find the pivot row
        for i in range(col, rows):
            if echelon_form[i, col] != 0:
                pivot_row = i
                break

        # If no nonzero pivot is found, move to the next column
        if pivot_row == -1:
            continue

        # Swap rows to move the pivot to the current row
        echelon_form.swap_rows(col, pivot_row)

        # Normalize the pivot row
        echelon_form[col, :] /= echelon_form[col, col]

        # Eliminate entries below the pivot
        for i in range(col + 1, rows):
            echelon_form[i, :] -= echelon_form[i, col] * echelon_form[col, :]

    return echelon_form

# Example Usage
# Create a matrix over the Rational Field
A = Matrix(QQ, [
    [1,  1, 1, 3],
    [2,  3, 7, 0],
    [1,  3, -2, 17]
])

# Apply Gaussian Elimination
echelon_form_A = gaussian_elimination(A)

# Display the result
print("Original Matrix A:")
print(A)
print("\nEchelon Form of A :")
print(echelon_form_A)

