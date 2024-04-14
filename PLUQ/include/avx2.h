#ifndef AVX2_H
#define AVX2_H

#include <stdio.h>
#include <immintrin.h>
#include "matrix.h"
#include "finite_field.h"

__m128i sub_avx2(__m128i a, __m128i b, __m128i vp);
long long reduce_vector2_int(__m256i input);
int scalar_product_avx2(const int *a, const int *b, int len, int p);
int scalar_product_scalar(int* v, int* w, int len, int p);
__m256d mul_mod_p(__m256d x, __m256d y, __m256d u, __m256d p);
void vector_scalar_mult_avx2(int *v, int size, int p, int c);
void vector_scalar_mult_scalar(int *v, int size, int p, int c);
void rows_elimination_avx2(int *A_data,int n, int matrixRank, int c,int p, __m256d vp, __m256d vu ,__m128i vp_128 , int k);
int scalar_product1_avx2(int* A_data, int matrixRank, int j, int n, int p);
void vectorMatrixMultiplication_avx2(int* A_data, int n, int matrixRank, int* result, int p);
__m256i sub256_avx2(__m256i a, __m256i b, __m256i vp);
void updateMatrix1_avx2(int* A_data, int matrixRank, int* tmp, int n,  int p);
int scalar_product2_avx2(int* A_data, int matrixRank, int n,int j, int pivot, int p);
void updateMatrix2_avx2(int* A_data, int matrixRank, int* tmp, int n, int m, int nullity, int pivot, int invpivot, int p);
#endif