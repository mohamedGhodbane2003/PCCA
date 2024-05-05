#ifndef UTILITIES_H
#define UTILITIES_H

#include "PLUQ.h"
#include "permutations.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

int* createRange(int n);
void printArray(int* array, int size);
bool checkTriL(Matrix L);
bool checkTriU(Matrix U, int rank);
void checkManyPLUQ(int p, int max_iter, int algo);
void checkOnePLUQ(int p, int m, int n, bool print, int algo);

#endif
