# Gaussian Elimination High-Performance Implementation

This program implements the high-performance PLUQ algorithm using AVX2 and compares it to the FLINT library implementation.

## Compilation

To compile, run `make`. To clean binary files, run `make clean`.

## Testing

To run tests:
- Execute `./bin/test`. This will show you choices to test this implementation:

    ***Algorithms***
    1) PLUQ
    2) PLUQ AVX2
    3) CROUT
    4) CROUT AVX2

- Run `./bin/test <choice>` to test the algorithm using all bit lengths (2 to 30).
- Run `./bin/test <choice> <bitlength> <m> <n> <0 or 1 (print)>` to test a specific case.

## Benchmarks

To run benchmarks:
- Execute `./bin/bench`. This will show you choices for the bench:

    ***Algorithms***
    1) PLUQ
    2) PLUQ AVX2
    3) CROUT
    4) CROUT AVX2
    5) FLINT

- Run `./bin/bench <choice> <size>` to benchmark the algorithm using all bit lengths (2 - 30) for a given size.
- Run `./bin/bench <choice> <size> <bitlength>` to benchmark the algorithm using a specific bit length for a given size.
- Run `./bin/bench 6` to benchmark all algorithms using different sizes (100, 300, 500, 1000, 1200) and different bit lengths (5, 12, 24, 29, 30), and store the result in bench.txt.
