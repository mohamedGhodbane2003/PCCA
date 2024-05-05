#include "avx2.h"

// Perform element-wise subtraction of two 128-bit integer vectors using AVX2 instructions
inline __m128i sub_avx2(__m128i a, __m128i b, __m128i vp) {
    // Subtract the corresponding elements of vectors a and b
    __m128i result = _mm_sub_epi32(a, b);
    // Create a mask to identify negative elements in the result vector
    __m128i mask = _mm_cmplt_epi32(result, _mm_setzero_si128());
    // Adjust the negative elements by adding the constant value p
    __m128i adjusted_result = _mm_add_epi32(result, vp);
    // Blend the adjusted result with the original result based on the mask
    result = _mm_blendv_epi8(result, adjusted_result, mask);
    return result;
}

// Perform element-wise subtraction of two 256-bit integer vectors using AVX2 instructions
// The same as '__m128i sub_avx2' with 256-bit
inline __m256i sub256_avx2(__m256i a, __m256i b, __m256i vp) {
    __m256i result = _mm256_sub_epi32(a, b);
    __m256i mask = _mm256_cmpgt_epi32(_mm256_setzero_si256(), result);
    __m256i adjusted_result = _mm256_add_epi32(result, vp);
    result = _mm256_blendv_epi8(result, adjusted_result, mask);
    return result;
}

// Function to reduce 4 elements of an AVX register (int64_t) to a single long long
long long reduce_vector2_int(__m256i input) {
      // Cast the AVX register to an array of long long and return the sum of the first 4 elements
    return ((long long*)&input)[0] + ((long long*)&input)[1] + ((long long*)&input)[2] + ((long long*)&input)[3];
}

// Function to compute the scalar product using AVX2 intrinsics
int scalar_product_avx2(int* A_data, int matrixRank, int k, int j, int n, int p) {
    // Initialize a vector to hold the sum of partial products (set to all zeros)
    __m256i sum_vec = _mm256_setzero_si256();
    int i;
    // Loop through elements in chunks of 4 (vector width)
    for(i = 0; i+3 < matrixRank; i += 4) {
        // Load 4 elements each from vectors 'A_data[k*n]' and 'A_data[i*n]' into AVX registers x and y
        __m256i x = _mm256_set_epi64x(A_data[k * n + (i+3)], A_data[k * n + (i+2)], A_data[k * n + (i+1)], A_data[k * n + i]);
        __m256i y = _mm256_set_epi64x(A_data[(i+3) * n + j], A_data[(i+2) * n + j], A_data[(i+1) * n + j], A_data[i * n + j]);
        // Perform element-wise multiplication between x and y and store the result in 'result'
        __m256i result =   _mm256_mul_epi32(x, y);
        // Add the multiplied elements in result to the running sum in sum_vec
        sum_vec = _mm256_add_epi64(sum_vec, result);
    }

    long long final = 0;
    // Handle the remaining elements (less than 4) using a scalar loop
    for(; i < matrixRank; ++i)
        final += (long long)A_data[k * n + i] * A_data[i * n + j];
    // Reduce the AVX register sum_vec to a single long long modulo p using the reduce_vector2_int function
    int res = (reduce_vector2_int(sum_vec) + final) % p;

    return res;
}

// Function to perform multiplication with modulo for double-precision floating-point vectors using AVX intrinsics
inline __m256d mul_mod_p(__m256d x, __m256d y, __m256d u, __m256d p) {
    __m256d h = _mm256_mul_pd(x, y); // Calculate high and low parts of the product (x * y)
    __m256d l = _mm256_fmsub_pd(x, y, h);
    __m256d b = _mm256_mul_pd(h, u);// Calculate intermediate value h/p
    __m256d c = _mm256_floor_pd(b); // integer part of b
    __m256d d = _mm256_fnmadd_pd(c, p, h);  // d = (c * p) - h (fused multiply-add)
    __m256d e = _mm256_add_pd(d, l);
    __m256d t = _mm256_sub_pd(e, p);
    // Conditional selection based on the sign of (e - p)
    e = _mm256_blendv_pd(t, e, t); // e = if (t < 0) then e else t
    // Handle potential overflow
    t = _mm256_add_pd(e, p);
    return _mm256_blendv_pd(e, t, e); // return e if (e < p) else t (conditional selection)
}

// Function to perform row elimination using AVX2 intrinsics
void rows_elimination_avx2(int *A_data,int n, int matrixRank, int c,int p , __m256d vp, __m256d vu ,__m128i vp_128, int k) {
    __m256d vc = _mm256_set1_pd(c); // - vc: constant vector holding the value 'c' (converted to double-precision)
    __m256d tmp; // - tmp: temporary variable to hold intermediate double-precision result
    int i;
    // Loop through elements in chunks of 4 (AVX register width)
    for (i = matrixRank + 1; i + 3 < n; i += 4) {
        // Load 4 integer elements from A_data into 128-bit AVX registers
        __m128i v1 = _mm_loadu_si128((__m128i *)&A_data[matrixRank * n + i]);
        __m128i v2 = _mm_loadu_si128((__m128i *)&A_data[k * n + i]);
        // Convert the loaded 128-bit integer vector (v1) to a 256-bit double-precision vector (vDouble)
        __m256d vDouble = _mm256_cvtepi32_pd(v1);
        // Calculate intermediate value using the mul_mod_p function for element-wise multiplication with modulo
        tmp = mul_mod_p(vc, vDouble, vu, vp);
        // Convert the double-precision result (tmp) back to a 128-bit integer vector
        __m128i resultInt = _mm256_cvttpd_epi32(tmp);
         // Perform vectorized subtraction with modulo
        __m128i result = sub_avx2(v2, resultInt, vp_128);
        // Store the updated result back into row 'k' of A_data
        _mm_storeu_si128((__m128i *)&A_data[k * n + i], result);
    }

    // This loop handles elements that don't fit into chunks of 4
     for (; i < n; ++i)
        A_data[k * n + i] = sub(A_data[k * n + i], mult(A_data[k * n + matrixRank], A_data[matrixRank * n + i], p), p);

}

// Function to perform vector-matrix multiplication using AVX2 intrinsics
void vectorMatrixMultiplication_avx2(int* A_data, int n, int matrixRank, int* result, int p){
        for(int i = matrixRank; i < n; i++){
            if(matrixRank > 3) // Check if the matrix rank (number of columns) is greater than 3 (vector width)
                // call scalar_product_avx2 function to compute the element (i-matrixRank)
                result[i-matrixRank] = scalar_product_avx2(A_data, matrixRank, matrixRank, i, n, p);
            else{// If matrix rank is small (less than or equal to 3), use a scalar loop
                result[i-matrixRank] = 0;
                for(int j = 0; j < matrixRank; j++)
                    result[i-matrixRank] = add(result[i-matrixRank], mult(A_data[matrixRank * n + j], A_data[j * n + i], p), p);
            } 
        }
}

// Function to update a matrix used in the Crout method
void updateMatrix1_avx2(int* A_data, int matrixRank, int* tmp, int n,  int p){
    __m256i vp = _mm256_set1_epi32(p); // - vp: a 256-bit AVX register set to all elements equal to 'p'
    int j;
    // This loop processes elements in groups of 8 to leverage AVX vectorized operations
    for (j = matrixRank; j + 7 < n; j += 8) {
         __m256i vA_data = _mm256_loadu_si256((__m256i*)(&A_data[matrixRank * n + j]));
         __m256i vtmp = _mm256_loadu_si256((__m256i*)(&tmp[-matrixRank + j]));
        // Perform element-wise subtraction between vA_data and vtmp with modulo 'p' using the sub256_avx2 function
         __m256i result = sub256_avx2(vA_data, vtmp, vp);
        // Store the updated result back into row 'matrixRank' of A_data (starting at index j)
         _mm256_storeu_si256((__m256i*)(&A_data[matrixRank * n + j]), result);
    }
    // This loop handles elements that don't fit into chunks of 8
    for(; j < n; j++)
        A_data[matrixRank * n + j] = sub(A_data[matrixRank * n + j], tmp[j-matrixRank], p);
}

// Function to update a matrix used in the Crout method
void updateMatrix2_avx2(int* A_data, int matrixRank, int* tmp, int n, int m, int nullity, int pivot, int invpivot, int p){
    // - vp: a 256-bit AVX register set to all elements equal to 'p' 
    __m256d vp = _mm256_set1_pd(p);
    // - vp_128: a 128-bit AVX register set to all elements equal to 'p'
    __m128i vp_128 = _mm_set1_epi32(p);
    // - vu: a 256-bit AVX register set to all elements equal to '1.0 / p'
    __m256d vu = _mm256_set1_pd(1.0 / p);
    // - vinvpivot: a 256-bit AVX register set to all elements equal to 'invpivot'
    __m256d vinvpivot = _mm256_set1_pd(invpivot);
    int tmp1_arr[4];
    int i;
    // This loop processes elements in groups of 4
    for(i = matrixRank + 1; i + 3 < m - nullity; i += 4){
        __m128i v1 = _mm_set_epi32(A_data[(i+3) * n + pivot], A_data[(i+2) * n + pivot], A_data[(i+1) * n + pivot], A_data[i * n + pivot]);
        __m128i v2 = _mm_loadu_si128((__m128i *)&tmp[i - matrixRank - 1]);
        // Perform element-wise subtraction between v1 and v2 with modulo 'p'
        __m128i tmp1 = sub_avx2(v1, v2, vp_128);
        // Convert tmp1 (128-bit integer) to a 256-bit double-precision vector for calculations involving mul_mod_p
        __m256d resultDouble = _mm256_cvtepi32_pd(tmp1);
        // Calculate an intermediate value using the mul_mod_p function
        resultDouble = mul_mod_p(resultDouble, vinvpivot, vu, vp);
        // Convert the updated double-precision result back to a 128-bit integer vector
        __m128i tmp2 = _mm256_cvttpd_epi32(resultDouble);
        // Store the temporary results from tmp2 into tmp1_arr
        _mm_storeu_si128((__m128i*)tmp1_arr, tmp2);
        // Update the corresponding elements in A_data (pivot column) using the values from tmp1_arr
        for(int j = 0; j < 4; j++)
            A_data[(i+j) * n + pivot] = tmp1_arr[j];
    }
    // This loop handles elements that don't fit into chunks of 4
    for(; i < m - nullity; i++){
        A_data[i * n + pivot] = sub(A_data[i * n + pivot], tmp[i - matrixRank - 1], p);
        A_data[i * n + pivot] = mult(A_data[i * n + pivot], invpivot, p);
    }
}

