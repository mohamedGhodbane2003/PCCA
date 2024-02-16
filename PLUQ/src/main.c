#include <stdio.h>
#include "utilities.h"

/** TODO: 
 *
 * Matrices:
 *
 * - no need for pointers to Matrix and malloc-ing Matrix, one should use
 *   Matrix structs directly for more efficiency (a Matrix is only two int's and a pointer)
 *
 * - we will target matrices of small dimensions (up to 100 or 200 rows/columns),
 *   in this case we will have better efficiency with a single contiguous storage space for all coefficients
 *   -> one malloc of an array of m*n entries
 *   -> entry (i,j) is stored at i*n + j in the array
 *
 *
 *
 *
 * Starting work on AVX2: use basic AVX intrinsics to write
 * . a function which takes two vectors v and w of the same length, filled with int32,
 * and a prime p, and computes the scalar product sum(v[i] w[i], i=0...len-1)
 * . a function which takes an int32 c, a vector v filled with int32, a prime p,
 * and computes the vector [c*v[0]  c*v[1]  ...  c*v[len-1]]
 *
 * Use all intrinsics you want from SSE and AVX:
 * https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#techs=SSE_ALL,AVX_ALL
 * (but not avx512 for the moment; or make one version with it and one version without it)
 *
 * Then compare the speed between the avx version and the usual version of the same function
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

