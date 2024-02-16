#include <stdio.h>
#include <time.h>
#include "utilities.h"

int main(int argc, char *argv[]) {
    // list of primes, primes[k] has bitlength k+2
    int primes[29] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457};

    if (argc != 4) {
        printf("Usage: %s <bitlength> <rows> <cols>\n", argv[0]);
        return 1;
    }

    int m = atoi(argv[2]);
    int n = atoi(argv[3]);
    int bitlength = atoi(argv[1]);

    if (bitlength < 2 || bitlength > 30) {
        printf("Invalid bitlength. Bitlength must be between 2 and 30.\n");
        return 1;
    }


    int p = primes[bitlength-2];

    // Create a matrix of specified dimensions
    Matrix *A = randomMatrix(m, n, p);

    int* P = NULL;
    int* Q = NULL;
    Matrix* LU = NULL;
    int rank = 0;
    
    // Benchmark PLUQ computation
    clock_t start = clock();
    PLUQ(A, &P, &LU, &Q, &rank, p);
    clock_t end = clock();

    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken to compute PLUQ: %f seconds\n", time_taken);

    // Cleanup
    destroyMatrix(A);

    return 0;
}
