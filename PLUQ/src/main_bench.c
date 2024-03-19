#include <stdio.h>
#include <time.h>
#include "utilities.h"
//#include <flint/fq_mat.h>
//#include <flint/fq.h>
//#include <flint/fmpz.h>
#include <flint/nmod_mat.h>
#include <flint/flint.h>


int main(int argc, char *argv[]) {
    int primes[29] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457};
    
    if (argc < 3 || (argv[1][0] != 'o' && argv[1][0] != 'f')) {
        printf("Usage: ./bin/bench <o/f> size\n");
        printf("  o: Benchmark AVX2 implementation against original implementation\n");
        printf("  f: Benchmark AVX2 implementation against FLINT implementation\n");
        return 1;
    }

    if (argc == 3 && (argv[1][0] == 'o' || argv[1][0] == 'f')) {
        int n = atoi(argv[2]);
        if(argv[1][0] == 'o'){
            printf("BITS\t\t\tORIGINAL PLUQ\t\t\tVECTORIZED PLUQ\t\t\tSPEEDUP\n\n");
            for (long k = 0; k < 29; k++) {
                double tt = 0.0;
                long nb_iter = 0;
                while (tt < 1.) {
                    Matrix A = randomMatrix(n, n, primes[k]);
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
                    Matrix A = randomMatrix(n, n, primes[k]);
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
        }else {
            printf("BITS\t\t\tFLINT PLUQ\t\t\tVECTORIZED PLUQ\t\t\tSPEEDUP\n\n");
                for (long k = 0; k < 29; k++) {
                    double tt = 0.0;
                    long nb_iter = 0;
                    while (tt < 1.) {
                        // VERSION FQ_MAT
                        //flint_rand_t state;
                        //fmpz_t p;
                        //fq_ctx_t ctx;
                        //fq_mat_t mat;
                        //slong rank;
                        //slong P[n];
                        //flint_randinit(state);
                        //flint_randseed(state, time(NULL), time(NULL));
                        //fmpz_init(p);
                        //fmpz_set_si(p, primes[k]);
                        //fq_ctx_init(ctx, p, 1, "x");
                        //fq_mat_init(mat, n, n, ctx);
                        //fq_mat_randtest(mat, state, ctx);
                        //clock_t start = clock();
                        //rank = fq_mat_lu(P , mat, rank, ctx);
                        //clock_t end = clock();
                        //tt += ((double)(end - start)) / CLOCKS_PER_SEC;

                        //fq_mat_clear(mat, ctx);
                        //fmpz_clear(p);
                        //flint_randclear(state);
                        //fq_ctx_clear(ctx);
                        //nb_iter += 1;

                        // VERSION NMOD_MAT
                        flint_rand_t state;
                        slong rank;
                        slong P[n];
                        flint_randinit(state);
                        flint_randseed(state, time(NULL), time(NULL));
                        nmod_mat_t mat;
                        nmod_mat_init(mat, n, n, primes[k]);
                        //nmod_mat_randtest(mat, state); --> the obtained matrix is too irregular
                        for (ulong i = 0; i < n; i++)
                            for (ulong j = 0; j < n; j++)
                                nmod_mat_entry(mat, i, j) = n_randint(state, primes[k]);
                        clock_t start = clock();
                        rank = nmod_mat_lu(P , mat, 0);  // 0 value for rank_check, otherwise there may be an early_exit
                        clock_t end = clock();
                        tt += ((double)(end - start)) / CLOCKS_PER_SEC;

                        nmod_mat_clear(mat);
                        flint_randclear(state);
                        nb_iter += 1;
                    }
            tt = tt / nb_iter;
            printf("%d\t\t\t%.12f\t\t\t",k + 2, tt);
            double tt_flint = tt;
            tt = 0.0;
            nb_iter = 0;
            while (tt < 1.) {
                Matrix A = randomMatrix(n, n, primes[k]);
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
            printf("%.12f\t\t\t%.4f\n\n", tt, tt_flint/tt);
            }
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


    double tt = 0.0;
    long nb_iter = 0;
    while (tt < 1.)
    {
        
        Matrix A = randomMatrix(m, n, p);

        int* P = NULL;
        int* Q = NULL;
        int rank = 0;

        
        clock_t start = clock();
        pluq_inplace_avx2(&A, &P, &Q, &rank, p);
        clock_t end = clock();
        tt += ((double)(end - start)) / CLOCKS_PER_SEC;

        
        free(A.data);

        nb_iter += 1;
    }

    tt = tt / nb_iter;
    printf("Time taken to compute PLUQ: %f seconds\n", tt);

    return 0;
}
