#include "avx2.h"

inline __m128i sub_avx2(__m128i a, __m128i b, int p) {
    __m128i vp = _mm_set1_epi32(p);
    __m128i result = _mm_sub_epi32(a, b);
    __m128i mask = _mm_cmplt_epi32(result, _mm_setzero_si128());
    __m128i adjusted_result = _mm_add_epi32(result, vp);
    result = _mm_blendv_epi8(result, adjusted_result, mask);
    return result;
}

long long reduce_vector2_int(__m256i input) {
    return ((long long*)&input)[0] + ((long long*)&input)[1] + ((long long*)&input)[2] + ((long long*)&input)[3];
}

int scalar_product_avx2(const int *a, const int *b, int len, int p) {
    __m256i sum_vec = _mm256_setzero_si256();

    for(int i = 0; i < len; i += 4) {
        __m256i x = _mm256_set_epi64x(a[i+3], a[i+2], a[i+1], a[i]);
        __m256i y = _mm256_set_epi64x(b[i+3], b[i+2], b[i+1], b[i]);
        __m256i result =   _mm256_mul_epi32(x, y);
        sum_vec = _mm256_add_epi64(sum_vec, result);
    }

    long long final = 0;
    for(int i = len - len % 4; i < len; ++i)
        final += (long long)a[i] * b[i];
    
    int res = (reduce_vector2_int(sum_vec) + final) % p;

    return res;
}


int scalar_product_scalar(int* v, int* w, int len, int p) {
    long long sum = 0;

    for (int i = 0; i < len; ++i) {
        sum += (long long)v[i] * w[i];
    }

    return sum % p;
}


__m256d mul_mod_p(__m256d x, __m256d y, __m256d u, __m256d p) {
    __m256d h = _mm256_mul_pd(x, y);
    __m256d l = _mm256_fmsub_pd(x, y, h);
    __m256d b = _mm256_mul_pd(h, u);
    __m256d c = _mm256_floor_pd(b);
    __m256d d = _mm256_fnmadd_pd(c, p, h);
    __m256d e = _mm256_add_pd(d, l);
    __m256d t = _mm256_sub_pd(e, p);
    e = _mm256_blendv_pd(t, e, t);
    t = _mm256_add_pd(e, p);
    return _mm256_blendv_pd(e, t, e);
}

void vector_scalar_mult_avx2(int *v, int size, int p, int c) {
    __m256d vc = _mm256_set1_pd(c);
    __m256d vp = _mm256_set1_pd(p);
    __m256d vu = _mm256_set1_pd(1.0 / p);
    __m256d result;
    for (int i = 0; i < size; i += 4) {
        // Load 4 consecutive 32-bit integers
        __m128i vInt = _mm_loadu_si128((__m128i *)&v[i]);

        // Convert 32-bit integers to double precision floating point numbers
        __m256d vDouble = _mm256_cvtepi32_pd(vInt);

        // Perform modular multiplication
        result = mul_mod_p(vc, vDouble, vu, vp);

        // Convert the result back to integers
        __m128i resultInt = _mm256_cvttpd_epi32(result);

        // Store the result back to memory
        _mm_storeu_si128((__m128i *)&v[i], resultInt);
    }
    long long tmp;
    for(int i = size - size % 4; i < size; ++i){
        tmp = (long long)v[i] * c;
        v[i] = tmp % p;
    }
}

void vector_scalar_mult_scalar(int *v, int size, int p, int c){
    long long tmp;
    for(int i = 0; i < size; i++){
        tmp = (long long) v[i] * c;
        v[i] = tmp % p;
    }
}

void print_m128i(__m128i vec) {
    int values[4];
    _mm_storeu_si128((__m128i*)values, vec);
    printf("[%d, %d, %d, %d]\n", values[0], values[1], values[2], values[3]);
}

void print_m256d(__m256d vec) {
    double values[4];
    _mm256_storeu_pd(values, vec);
    printf("[%f, %f, %f, %f]\n", values[0], values[1], values[2], values[3]);
}

void rows_elimination_avx2(int *A_data,int n, int matrixRank, int c, int p, int k) {
    __m256d vc = _mm256_set1_pd(c);
    __m256d vp = _mm256_set1_pd(p);
    __m256d vu = _mm256_set1_pd(1.0 / p);
    __m256d tmp;
    int i;
    for (i = matrixRank + 1; i + 3 < n; i += 4) {
        __m128i v1 = _mm_loadu_si128((__m128i *)&A_data[matrixRank * n + i]);
        __m128i v2 = _mm_loadu_si128((__m128i *)&A_data[k * n + i]);
        __m256d vDouble = _mm256_cvtepi32_pd(v1);
        tmp = mul_mod_p(vc, vDouble, vu, vp);
        __m128i resultInt = _mm256_cvttpd_epi32(tmp);
        __m128i result = sub_avx2(v2, resultInt, p);
        _mm_storeu_si128((__m128i *)&A_data[k * n + i], result);
    }

     for (; i < n; ++i)
        A_data[k * n + i] = sub(A_data[k * n + i], mult(A_data[k * n + matrixRank], A_data[matrixRank * n + i], p), p);

}

