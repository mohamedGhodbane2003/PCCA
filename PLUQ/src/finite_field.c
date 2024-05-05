#include "finite_field.h"

// Adds two 32-bit integers 'a' and 'b' and returns the result modulo 'p'
inline int add(int a, int b, int p) {
    int r = a + b;
    return r >= p ?  r - p : r;
}

// Subtracts 32-bit integer 'b' from 'a' and returns the result modulo 'p'
inline int sub(int a, int b, int p) {
    int r = a - b;
    return r < 0 ?  r + p : r;
}

// multiplies two 32-bit integers 'a' and 'b' and returns the result modulo 'p'
inline int mult(int a, int b, int p) {
    // Use 64-bit integer type to avoid overflow
    long long r = (long long)a * b;
    return r % p;
}

// Computes the modular multiplicative inverse of 'a' modulo 'p'
inline int inverse(int a, int p) {
    int m0 = p, t, q;
    int x0 = 0, x1 = 1;
    if (p == 1)
       return 0;
    // extended Euclid Algorithm
    while (a > 1) {
        // q is quotient
        q = a / p;
        t = p;
        // p is remainder now, process same as Euclid's algo
        p = a % p, a = t;
        t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    // Make it positive
    if (x1 < 0)
       x1 += m0;
    return x1;
}

