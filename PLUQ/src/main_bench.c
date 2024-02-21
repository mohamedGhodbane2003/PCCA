#include <stdio.h>
#include <time.h>
#include "utilities.h"

int main(int argc, char *argv[]) {
    // list of primes, primes[k] has bitlength k+2
    int primes[29] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457};

    if (argc != 4) {
        printf("BITS\t\t\tORIGINAL PLUQ\t\t\tVECTORIZED PLUQ\t\t\tSPEEDUP\n\n");
        for (long k = 0; k < 29; k++) {
            double tt = 0.0;
            long nb_iter = 0;
            while (tt < 1.) {
                Matrix A = randomMatrix(200, 200, primes[k]);
                int* P = NULL;
                int* Q = NULL;
                int rank = 0;
                clock_t start = clock();
                pluq_inplace(&A, &P, &Q, &rank, primes[k]);
                clock_t end = clock();
                tt += ((double)(end - start)) / CLOCKS_PER_SEC;
                
                nb_iter += 1;
            }

            tt = tt / nb_iter;
            printf("%d\t\t\t%.12f\t\t\t",k + 2, tt);
            double tt_original = tt;
            tt = 0.0;
            nb_iter = 0;
            while (tt < 1.) {
                Matrix A = randomMatrix(200, 200, primes[k]);
                int* P = NULL;
                int* Q = NULL;
                int rank = 0;
                clock_t start = clock();
                pluq_inplace_avx2(&A, &P, &Q, &rank, primes[k]);
                clock_t end = clock();
                tt += ((double)(end - start)) / CLOCKS_PER_SEC;
                
                nb_iter += 1;
            }

            tt = tt / nb_iter;
            printf("%.12f\t\t\t%.4f\n\n", tt, tt_original/tt);
        }
        return 0;
    }

    int m = atoi(argv[2]);
    int n = atoi(argv[3]);
    int bitlength = atoi(argv[1]);

    if (bitlength < 2 || bitlength > 30) {
        printf("Invalid bitlength. Bitlength must be between 2 and 30.\n");
        return 1;
    }


    int p = primes[bitlength-2];

    // repeat computation several times if it is < 1 second
    double tt = 0.0;
    long nb_iter = 0;
    while (tt < 1.)
    {
        // Create a matrix of specified dimensions
        Matrix A = randomMatrix(m, n, p);

        int* P = NULL;
        int* Q = NULL;
        Matrix LU;
        int rank = 0;

        // Benchmark PLUQ computation
        clock_t start = clock();
        PLUQ(A, &P, &LU, &Q, &rank, p);
        clock_t end = clock();
        tt += ((double)(end - start)) / CLOCKS_PER_SEC;

        // Cleanup
        free(A.data);

        nb_iter += 1;
    }

    tt = tt / nb_iter;
    printf("Time taken to compute PLUQ: %f seconds\n", tt);

    return 0;
}
