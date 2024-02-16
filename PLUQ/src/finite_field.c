#include "finite_field.h"

/*
 * Function: add
 * --------------
 * Performs addition in the finite field.
 * 
 * Parameters:
 *   - a: First operand.
 *   - b: Second operand.
 *   - p: Prime number used for modulo operation.
 * 
 * Returns:
 *   The result of the addition operation modulo p.
 */

inline int add(int a, int b, int p) {
    int r = a + b;
    // Reduce the result modulo p
    return r >= p ?  r - p : r;
}

/*
 * Function: sub
 * --------------
 * Performs subtraction in the finite field.
 * 
 * Parameters:
 *   - a: First operand.
 *   - b: Second operand.
 *   - p: Prime number used for modulo operation.
 * 
 * Returns:
 *   The result of the subtraction operation modulo p.
 */

inline int sub(int a, int b, int p) {
    int r = a - b;
    // Reduce the result modulo p
    return r < 0 ?  r + p : r;
}

/*
 * Function: mult
 * ---------------
 * Performs multiplication in the finite field.
 * 
 * Parameters:
 *   - a: First operand.
 *   - b: Second operand.
 *   - p: Prime number used for modulo operation.
 * 
 * Returns:
 *   The result of the multiplication operation modulo p.
 */

inline int mult(int a, int b, int p) {
    // Use 64-bit integer type to avoid overflow
    long long r = (long long)a * b;
    // Reduce the result modulo p
    return r % p;
}

/*
 * Function: inverse
 * ------------------
 * Calculates the multiplicative inverse of a modulo p.
 * 
 * Parameters:
 *   - a: The number for which the inverse is to be found.
 *   - p: Prime number.
 * 
 * Returns:
 *   The multiplicative inverse of a modulo p.
 */

int inverse(int a, int p) {
    // Using extended Euclidean algorithm
    int m0 = p, t, q;
    int x0 = 0, x1 = 1;

    if (p == 1)
       return 0;

    // Apply extended Euclid Algorithm
    while (a > 1) {
        // q is quotient
        q = a / p;
        t = p;

        // p is remainder now, process same as
        // Euclid's algo
        p = a % p, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }

    // Make x1 positive
    if (x1 < 0)
       x1 += m0;

    return x1;
}

