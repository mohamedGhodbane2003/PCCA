#include <stdio.h>
#include <time.h>
#include "utilities.h"
#include <flint/nmod_mat.h>
#include <flint/flint.h>


int main(int argc, char *argv[]) {
    int primes[29] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289,
                      24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 
                      12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457};
    
    if (argc == 1) {
        printf("***Algorithms***\n");
        printf("1) PLUQ\n");
        printf("2) PLUQ AVX2\n");
        printf("3) CROUT\n");
        printf("4) CROUT AVX2\n");
        printf("5) FLINT\n");
        printf("./bin/bench <choice> <size> : bench the algorithm using all bitlengths (2 - 30) \n");
        printf("./bin/bench <choice> <size> <bitlength> : bench the algorithm using <bitlength> \n");
        printf("./bin/bench 6 : bench all the algorithm using diffrent sizes (100, 300, 500, 1000, 1200) and diffrent bitlengths (5, 12, 24, 29, 30) and store the result in bench.txt\n");
        return 1;
    }

  if (atoi(argv[1]) == 6 && argc == 2){
    int sizes[5] = {100, 300, 500, 1000, 1200};
    int bitlengths[5] = {5, 12, 24, 29, 30};
     FILE *file = fopen("bench.txt", "w"); 
        if (file == NULL) {
            perror("Error opening file");
            return 1;
        }
    for(int i = 4; i < 5; i++ ){
        fprintf(file, "Prime Bitlength\t\t\tSize\t\t\t");
        switch(i+1){
            case 1: fprintf(file, "PLUQ\n\n"); printf("bench PLUQ...\n"); break;
            case 2: fprintf(file, "PLUQ AVX2\n\n"); printf("bench PLUQ AVX2...\n"); break;
            case 3: fprintf(file, "Crout\n\n"); printf("bench Crout...\n"); break;
            case 4: fprintf(file, "Crout AVX2\n\n"); printf("bench Crout AVX2...\n"); break;
            case 5: fprintf(file, "FLINT\n\n"); printf("bench FLINT...\n"); break;
            default: break;
        }
        for(int j  = 0; j < 5; j++){
            for(int k = 0; k < 5; k++){
                printf("%d, %d\n", j, k);
               double tt = 0.0;
            long nb_iter = 0;
            while (tt < 1.) {
                // VERSION NMOD_MAT
                flint_rand_t state;
                slong rankF;
                slong PF[sizes[j]];
                flint_randinit(state);
                flint_randseed(state, time(NULL), time(NULL));
                nmod_mat_t mat;
                nmod_mat_init(mat, sizes[j], sizes[j], primes[bitlengths[k - 2]]);
                //nmod_mat_randtest(mat, state); --> the obtained matrix is too irregular
                for (ulong l = 0; l < sizes[j]; l++)
                    for (ulong m = 0; m < sizes[j]; m++)
                        nmod_mat_entry(mat, l, m) = n_randint(state, primes[bitlengths[k - 2]]);
                Matrix A = randomMatrix(sizes[j], sizes[j], primes[bitlengths[k - 2]]);
                int* P = NULL;
                int* Q = NULL;
                int rank = 0;
                clock_t start, end;
                switch(i+1){
                    case 1: start = clock();
                            pluq_inplace(&A, &P, &Q, &rank, primes[bitlengths[k - 2]]);
                            end = clock(); break;
                    case 2: start = clock();
                            pluq_inplace_avx2(&A, &P, &Q, &rank, primes[bitlengths[k - 2]]);
                            end = clock(); break;
                    case 3: start = clock();
                            pluq_crout(&A, &P, &Q, &rank, primes[bitlengths[k - 2]]);
                            end = clock(); break;
                    case 4: start = clock();
                            pluq_crout_avx2(&A, &P, &Q, &rank, primes[bitlengths[k - 2]]);
                            end = clock(); break;
                    case 5: start = clock();
                            rankF = nmod_mat_lu(PF , mat, 0);  // 0 value for rank_check, otherwise there may be an early_exit
                            end = clock(); break;
                    default: break;
                }
                tt += ((double)(end - start)) / CLOCKS_PER_SEC;
                nb_iter++;
                nmod_mat_clear(mat);
                flint_randclear(state);
            }
            tt = tt / nb_iter;
            fprintf(file, "%d\t\t\t\t%d\t\t\t%.12f\n", bitlengths[k], sizes[j], tt);
            }
        }
        fprintf(file, "\n\n\n\n");
    }
    fclose(file);
    printf("Results in file bench.txt\n");
    return 0;
  } 
  if (argc == 4){
        int choice = atoi(argv[1]);
        int n = atoi(argv[2]);
        int bitlength = atoi(argv[3]);
        printf("Prime Bitlength\t\t\t");
        switch(choice){
            case 1: printf("PLUQ\n\n"); break;
            case 2: printf("PLUQ AVX2\n\n"); break;
            case 3: printf("Crout\n\n"); break;
            case 4: printf("Crout AVX2\n\n"); break;
            case 5: printf("FLINT\n\n"); break;
            default: break;
        }
            double tt = 0.0;
            long nb_iter = 0;
            while (tt < 1.) {
                // VERSION NMOD_MAT
                flint_rand_t state;
                slong rankF;
                slong PF[n];
                flint_randinit(state);
                flint_randseed(state, time(NULL), time(NULL));
                nmod_mat_t mat;
                nmod_mat_init(mat, n, n, primes[bitlength - 2]);
                //nmod_mat_randtest(mat, state); --> the obtained matrix is too irregular
                for (ulong i = 0; i < n; i++)
                    for (ulong j = 0; j < n; j++)
                        nmod_mat_entry(mat, i, j) = n_randint(state, primes[bitlength - 2]);
                Matrix A = randomMatrix(n, n, primes[bitlength - 2]);
                int* P = NULL;
                int* Q = NULL;
                int rank = 0;
                clock_t start, end;
                switch(choice){
                    case 1: start = clock();
                            pluq_inplace(&A, &P, &Q, &rank, primes[bitlength - 2]);
                            end = clock(); break;
                    case 2: start = clock();
                            pluq_inplace_avx2(&A, &P, &Q, &rank, primes[bitlength - 2]);
                            end = clock(); break;
                    case 3: start = clock();
                            pluq_crout(&A, &P, &Q, &rank, primes[bitlength - 2]);
                            end = clock(); break;
                    case 4: start = clock();
                            pluq_crout_avx2(&A, &P, &Q, &rank, primes[bitlength - 2]);
                            end = clock(); break;
                    case 5: start = clock();
                            rankF = nmod_mat_lu(PF , mat, 0);  // 0 value for rank_check, otherwise there may be an early_exit
                            end = clock(); break;
                    default: break;
                }
                tt += ((double)(end - start)) / CLOCKS_PER_SEC;
                nb_iter++;
                nmod_mat_clear(mat);
                flint_randclear(state);
            }
            tt = tt / nb_iter;
            printf("%d\t\t\t%.12f\n", bitlength, tt);
        return 0;
  }

  if (argc == 3) {
        int choice = atoi(argv[1]);
        int n = atoi(argv[2]);
        printf("Prime Bitlength\t\t\t");
        switch(choice){
            case 1: printf("PLUQ\n\n"); break;
            case 2: printf("PLUQ AVX2\n\n"); break;
            case 3: printf("Crout\n\n"); break;
            case 4: printf("Crout AVX2\n\n"); break;
            case 5: printf("FLINT\n\n"); break;
            default: break;
        }
        for (long k = 0; k < 29; k++) {
            double tt = 0.0;
            long nb_iter = 0;
            while (tt < 1.) {
                // VERSION NMOD_MAT
                flint_rand_t state;
                slong rankF;
                slong PF[n];
                flint_randinit(state);
                flint_randseed(state, time(NULL), time(NULL));
                nmod_mat_t mat;
                nmod_mat_init(mat, n, n, primes[k]);
                //nmod_mat_randtest(mat, state); --> the obtained matrix is too irregular
                for (ulong i = 0; i < n; i++)
                    for (ulong j = 0; j < n; j++)
                        nmod_mat_entry(mat, i, j) = n_randint(state, primes[k]);
                Matrix A = randomMatrix(n, n, primes[k]);
                int* P = NULL;
                int* Q = NULL;
                int rank = 0;
                clock_t start, end;
                switch(choice){
                    case 1: start = clock();
                            pluq_inplace(&A, &P, &Q, &rank, primes[k]);
                            end = clock(); break;
                    case 2: start = clock();
                            pluq_inplace_avx2(&A, &P, &Q, &rank, primes[k]);
                            end = clock(); break;
                    case 3: start = clock();
                            pluq_crout(&A, &P, &Q, &rank, primes[k]);
                            end = clock(); break;
                    case 4: start = clock();
                            pluq_crout_avx2(&A, &P, &Q, &rank, primes[k]);
                            end = clock(); break;
                    case 5: start = clock();
                            rankF = nmod_mat_lu(PF , mat, 0);  // 0 value for rank_check, otherwise there may be an early_exit
                            end = clock(); break;
                    default: break;
                }
                tt += ((double)(end - start)) / CLOCKS_PER_SEC;
                nb_iter++;
                nmod_mat_clear(mat);
                flint_randclear(state);
            }
            tt = tt / nb_iter;
            printf("%d\t\t\t%.12f\n", k + 2, tt);
        }
        return 0;
    }
    return 0;
}
