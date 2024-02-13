#include <stdio.h>
#include "utilities.h"

/** TODO: 
 *
 *  \todo this does not seem to work for primes of 16 bits or more
 *  (overflow issue in mult or inverse in finitefield.c ?)
 *
 *  \todo for TEST:
 *  via argv / argc, propose the choice of a bitlength to the user for the prime number
 *     --> if choice given, take prime of that length in the list, and check
 *     --> if no choice given, check all bitlengths like is already done now
 *  ./bin/test 5 10 20     --> launches test on 10 x 20 matrix, with prime of bitlength 5
 *  ./bin/test             --> launches test on small matrices and many primes (like already done now)
 *
 *  \todo add benchmarks:
 *  the user gives number of rows, number of columns, and this prints the time for computing PLUQ
 *  ./bin/bench 5 10 20     --> launches benchmark on 10 x 20 matrix, with prime of bitlength 5
 *
 *  This could be done by two different executables from two files main_test.c and main_bench.c, managed by the Makefile (make test || make bench)
 */

/*
 * Note: Overflow issue was in multiplyMatrices function in file matrix.c has been resolved
 *       by casting intermediate results to 'long long'
 */

int main(int argc, char *argv[]) {
    // list of primes, primes[k] has bitlength k+2
    int primes[29] = {3, 7, 13, 29, 53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189, 805306457};

    // Check if a bit length, m and n are provided as command-line argument
    if (argc == 4) {

        int bitlength = atoi(argv[1]);
        int m = atoi(argv[2]);
        int n = atoi(argv[3]);
        
        if (bitlength < 2 || bitlength > 30) {
            printf("Invalid bitlength. Bitlength must be between 2 and 30.\n");
            return 1;
        }
        int p = primes[bitlength - 2];
        printf("Prime %d, %d bits\n", p, bitlength);
        checkOnePLUQ(p, m, n, false);
        printf("\n");
        return 0;
    }

    // If no bit length is provided, test for all primes in the list
    for (long k = 0; k < 29; k++) {
        printf("Prime %d, %ld bits\n", primes[k], k + 2);
        checkManyPLUQ(primes[k], 1000);
        printf("\n");
    }

    return 0;
}
/** computation of primes of different lengths
sage: for nbits in range(4,31):
....:     n = 2**(nbits-1) + 2**(nbits-2)
....:     print(f"{n.next_prime()}      \t// {nbits} bits")
....:
13              // 4 bits
29              // 5 bits
53              // 6 bits
97              // 7 bits
193             // 8 bits
389             // 9 bits
769             // 10 bits
1543            // 11 bits
3079            // 12 bits
6151            // 13 bits
12289           // 14 bits
24593           // 15 bits
49157           // 16 bits
98317           // 17 bits
196613          // 18 bits
393241          // 19 bits
786433          // 20 bits
1572869         // 21 bits
3145739         // 22 bits
6291469         // 23 bits
12582917        // 24 bits
25165843        // 25 bits
50331653        // 26 bits
100663319       // 27 bits
201326611       // 28 bits
402653189       // 29 bits
805306457       // 30 bits
*/

