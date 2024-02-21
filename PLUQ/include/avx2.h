#ifndef AVX2_H
#define AVX2_H

#include <stdio.h>
#include <immintrin.h>
#include "matrix.h"
#include "finite_field.h"

__m128i sub_avx2(__m128i a, __m128i b, int p);
long long reduce_vector2_int(__m256i input);
int scalar_product_avx2(const int *a, const int *b, int len, int p);
int scalar_product_scalar(int* v, int* w, int len, int p);
__m256d mul_mod_p(__m256d x, __m256d y, __m256d u, __m256d p);
void vector_scalar_mult_avx2(int *v, int size, int p, int c);
void vector_scalar_mult_scalar(int *v, int size, int p, int c);
void rows_elimination_avx2(int *A_data,int n, int matrixRank, int c, int p, int k);

#endif