#ifndef PLUQ_H
#define PLUQ_H

#include "matrix.h"
#include "utilities.h"
#include "permutations.h"
#include "finite_field.h"
#include "avx2.h"
#include <stdio.h>
#include <stdlib.h>

void PLUQ(Matrix A, int** P, Matrix* LU, int** Q, int* rank, int p);
void expand_PLUQ(Matrix LU, int rank, Matrix* L, Matrix* U);
void pluq_inplace(Matrix* A, int** P, int** Q, int* rank, int p);
void pluq_inplace_avx2(Matrix* A, int** P, int** Q, int* rank, int p);

#endif
