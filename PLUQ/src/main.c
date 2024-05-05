#include <stdio.h>
#include "utilities.h"

int main(int argc, char *argv[]) {
    // list of primes, primes[k] has bitlength k+2
    int primes[29] = {  3, 7, 13, 29, 53, 97, 193, 389, 769, 1543,
                        3079, 6151, 12289, 24593, 49157, 98317, 196613,
                        393241, 786433, 1572869, 3145739, 6291469, 12582917,
                        25165843, 50331653, 100663319, 201326611, 402653189, 805306457};

    // If no command-line arguments are provided, print the list of choices
    if (argc == 1) {
        printf("***Algorithms***\n");
        printf("1) PLUQ\n");
        printf("2) PLUQ AVX2\n");
        printf("3) CROUT\n");
        printf("4) CROUT AVX2\n");
        printf(".\\bin\\test <choice> : Test the algorithm using all bitlengths (2 to 30)\n.\\bin\\test <choice> <bitlength> <m> <n> <0 or 1 (print)> : to test specific case\n");
        return 0;
    }

    // test for all primes
    if(argc == 2){
        int choice = atoi(argv[1]);
        for (long k = 0; k < 29; k++) {
            printf("Prime %d, %ld bits\n", primes[k], k + 2);
            checkManyPLUQ(primes[k], 1000, choice);
            printf("\n");
        }
    }

    // Check if a bit length, m and n are provided as command-line argument
    if (argc == 6) {
        int choice = atoi(argv[1]);
        int bitlength = atoi(argv[2]);
        int m = atoi(argv[3]);
        int n = atoi(argv[4]);
        int print_flag = atoi(argv[5]);
        
        if (bitlength < 2 || bitlength > 30) {
            printf("Invalid bitlength. Bitlength must be between 2 and 30.\n");
            return 1;
        }
        int p = primes[bitlength - 2];
        printf("Prime %d, %d bits\n", p, bitlength);
        checkOnePLUQ(p, m, n, print_flag, choice);
        printf("\n");
        return 0;
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

