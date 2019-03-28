#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "getoptions.h"
#include "ktiming.h"

#ifndef RAND_MAX
#define RAND_MAX 32767
#endif

#ifndef TIMING_COUNT
#define TIMING_COUNT 0
#endif

#ifndef SERIAL
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#else
#define cilk_spawn
#define cilk_sync
#define __cilkrts_init()
#define __cilkrts_reset_timing()
#define __cilkrts_accum_timing()
#endif


#ifdef NO_PIN
#define __cilkrts_set_pinning_info(n)
#define __cilkrts_unset_pinning_info()
#define __cilkrts_pin_top_level_frame_at_socket(n)
#endif


/* The real numbers we are using --- either double or float */
typedef double REAL;

/* maximum tolerable relative error (for the checking routine) */
#define EPSILON (1.0E-6)
#define CACHE_LINE_SIZE 64

static const unsigned int POWER = 4;
static const unsigned int DAC_ARITH_BASECASE = (1 << POWER);  // 16x16
static const unsigned int MATMUL_THRESH = (1 << (POWER + 1)); // 64x64

int top_n;

/* n is the current matrix size of M, and
 * orig_n is the original matrix that M is part of
 */
#define Z_PARTITION(M, M1, M2, M3, M4, n, orig_n) \
  M1 = &M[0]; \
  M2 = &M[n >> 1]; \
  M3 = &M[(n * orig_n) >> 1]; \
  M4 = &M[((n * orig_n) + n) >> 1]; \

/*
 * Matrices are stored in row-major order; A is a pointer to
 * the first element of the matrix, and an is the number of elements
 * between two rows. This macro produces the element A[i,j]
 * given A, an, i and j
 */
#define ELEM(A, an, i, j) (A[(i) * (an) + (j)])

#define UNROLL(M, m) \
  REAL m##_0, m##_1, m##_2, m##_3, m##_4, m##_5, m##_6, m##_7; \
  m##_0 = *M; m##_1 = *(M+1); m##_2 = *(M+2); m##_3 = *(M+3); \
  m##_4 = *(M+4); m##_5 = *(M+5); m##_6 = *(M+6); m##_7 = *(M+7);

#define UNROLL_ADD_1(C, t) \
  (*(C))   += t##_0; \
  (*(C+1)) += t##_1; \
  (*(C+2)) += t##_2; \
  (*(C+3)) += t##_3; \
  (*(C+4)) += t##_4; \
  (*(C+5)) += t##_5; \
  (*(C+6)) += t##_6; \
  (*(C+7)) += t##_7;

#define UNROLL_ASSIGN(C, t) \
  (*(C))   = t##_0; \
  (*(C+1)) = t##_1; \
  (*(C+2)) = t##_2; \
  (*(C+3)) = t##_3; \
  (*(C+4)) = t##_4; \
  (*(C+5)) = t##_5; \
  (*(C+6)) = t##_6; \
  (*(C+7)) = t##_7;

#define UNROLL_ADD_3(C, u, v, w) \
  (*(C))   += u##_0 + v##_0 + w##_0; \
  (*(C+1)) += u##_1 + v##_1 + w##_1; \
  (*(C+2)) += u##_2 + v##_2 + w##_2; \
  (*(C+3)) += u##_3 + v##_3 + w##_3; \
  (*(C+4)) += u##_4 + v##_4 + w##_4; \
  (*(C+5)) += u##_5 + v##_5 + w##_5; \
  (*(C+6)) += u##_6 + v##_6 + w##_6; \
  (*(C+7)) += u##_7 + v##_7 + w##_7;

#define UNROLL_VAL_ADD_2(c, a, b) \
  c##_0 += a##_0 + b##_0; \
  c##_1 += a##_1 + b##_1; \
  c##_2 += a##_2 + b##_2; \
  c##_3 += a##_3 + b##_3; \
  c##_4 += a##_4 + b##_4; \
  c##_5 += a##_5 + b##_5; \
  c##_6 += a##_6 + b##_6; \
  c##_7 += a##_7 + b##_7; \

#define UNROLL_VAL_SUB_SELF(C, t) \
  (*(C))   = t##_0 - (*(C)); \
  (*(C+1)) = t##_1 - (*(C+1)); \
  (*(C+2)) = t##_2 - (*(C+2)); \
  (*(C+3)) = t##_3 - (*(C+3)); \
  (*(C+4)) = t##_4 - (*(C+4)); \
  (*(C+5)) = t##_5 - (*(C+5)); \
  (*(C+6)) = t##_6 - (*(C+6)); \
  (*(C+7)) = t##_7 - (*(C+7));

#define UNROLL_ADD_2_ASSIGN(C, a, b) \
  (*(C))   = a##_0 + b##_0; \
  (*(C+1)) = a##_1 + b##_1; \
  (*(C+2)) = a##_2 + b##_2; \
  (*(C+3)) = a##_3 + b##_3; \
  (*(C+4)) = a##_4 + b##_4; \
  (*(C+5)) = a##_5 + b##_5; \
  (*(C+6)) = a##_6 + b##_6; \
  (*(C+7)) = a##_7 + b##_7;

#define UNROLL_SUB_2_ASSIGN(C, a, b) \
  (*(C))   = a##_0 - b##_0; \
  (*(C+1)) = a##_1 - b##_1; \
  (*(C+2)) = a##_2 - b##_2; \
  (*(C+3)) = a##_3 - b##_3; \
  (*(C+4)) = a##_4 - b##_4; \
  (*(C+5)) = a##_5 - b##_5; \
  (*(C+6)) = a##_6 - b##_6; \
  (*(C+7)) = a##_7 - b##_7;

#define UNROLL_PRINT(C) \
  printf("%p\n", ((C))); \
  printf("%p\n", ((C+1))); \
  printf("%p\n", ((C+2))); \
  printf("%p\n", ((C+3))); \
  printf("%p\n", ((C+4))); \
  printf("%p\n", ((C+5))); \
  printf("%p\n", ((C+6))); \
  printf("%p\n", ((C+7)));


unsigned long rand_nxt = 0;

int cilk_rand(void) {
  int result;
  rand_nxt = rand_nxt * 1103515245 + 12345;
  result = (rand_nxt >> 16) % ((unsigned int) RAND_MAX + 1);
  return result;
}


/*
 * Set an n by n matrix A to random values.  The distance between
 * rows is an
 */
static void init_matrix(int n, REAL *A, int an) {

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      A[i * n + j] = i * n + j;//((double) cilk_rand()) / (double) RAND_MAX;
    }
  }
}

/*
 * Set an n by n matrix A to random values.  The distance between
 * rows is an
 */
static void print_matrix(int n, REAL *A, int orig_n) {

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      printf("%4.0f ", A[i * orig_n + j]);
    }
    printf("\n");
  }
}


/*
 * Naive sequential algorithm, for comparison purposes
 */
void matrixmul(REAL *C, REAL *A, REAL *B, int n, int orig_n) {

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      REAL s = (REAL)0;
      for(int k = 0; k < n; ++k) {
        s += A[i * orig_n + k] * B[k * orig_n + j];
      }
      C[i * orig_n + j] = s;
    }
  }
}

/**
 * Assumption: by the time we get here the matrices are within
 * the row-major order with row width exactly n
 * C: output, size nxn
 * A, B: input, size nxn
 * Same as mm_base except that it sums the original value of C into it
 **/
static void mm_additive_base(REAL *C, REAL *A, REAL *B, int n, int c_orig_n, int a_orig_n, int b_orig_n) {

  REAL *ptrToA = A;
  REAL *ptrToB = B;

  for(int row = 0; row < n; row++) { // going down the row
    for(int col = 0; col < n; col+=8) { // doing 8 columns at a time
      ptrToB = B;  // put ptrToB back to the beginning of B
      REAL valA = *ptrToA;
      ptrToA++; // advance A to the next element in the row block
      REAL c0 = *(C)   + valA * ptrToB[col];
      REAL c1 = *(C+1) + valA * ptrToB[col+1];
      REAL c2 = *(C+2) + valA * ptrToB[col+2];
      REAL c3 = *(C+3) + valA * ptrToB[col+3];
      REAL c4 = *(C+4) + valA * ptrToB[col+4];
      REAL c5 = *(C+5) + valA * ptrToB[col+5];
      REAL c6 = *(C+6) + valA * ptrToB[col+6];
      REAL c7 = *(C+7) + valA * ptrToB[col+7];

      // accumulate into c_i the product that needs to go into current row
      // block of C
      for(int j = 1; j < n; j++) {
        // each iter does a block of columns in the current row
        ptrToB += b_orig_n; // advance B to the next row
        valA = *ptrToA;
        ptrToA++; // advance A to the next element
        c0 += valA * ptrToB[col];
        c1 += valA * ptrToB[col+1];
        c2 += valA * ptrToB[col+2];
        c3 += valA * ptrToB[col+3];
        c4 += valA * ptrToB[col+4];
        c5 += valA * ptrToB[col+5];
        c6 += valA * ptrToB[col+6];
        c7 += valA * ptrToB[col+7];
      }
      ptrToA -= n; // move A back to beginning of the row
      *(C) = c0;
      *(C+1) = c1;
      *(C+2) = c2;
      *(C+3) = c3;
      *(C+4) = c4;
      *(C+5) = c5;
      *(C+6) = c6;
      *(C+7) = c7;
      //assert(C == &oldC[row*n + col]);
      C += 8;
    }
    //assert(C == &oldC[row*n + n]);
    C += (c_orig_n - n);
    ptrToA += a_orig_n; // advance A to the next row
  }
  return;
}

/**
 * Assumption: by the time we get here the matrices are within
 * the row-major order with row width exactly n, i.e., all nxn
 * elements are laid out in continuguous memory without breaks.
 * C: output, size nxn
 * A, B: input, size nxn
 *
 * This base is optimized to do a block of column
 * (8 contiguous elements) at a time.
 * It walks through each C[i,j] in row major order, but compute
 * C[i,j]--C[i+7,j] (inclusive) at once in each iteration of 2nd loop.
 * The way it does is as follows:
 * Each C[i,j] needs to sum together the product between every elements
 * in row i of A and col j of B.  Each inner most loop is effectively doing
 * the product necessary involving a single element in A and traverses
 * the entire row.
 * Each iteration in the 2nd loop traverses a given row of A over and over,
 * doing multiplication necessary involving the next row of column block.
 * does the nece
 **/
static void mm_base(REAL *C, REAL *A, REAL *B, int n, int c_orig_n, int a_orig_n, int b_orig_n) {

  REAL *ptrToA = A;
  REAL *ptrToB = B;
  int dummy = 0;
  for(int row = 0; row < n; row++) { // going down the row
    for(int col = 0; col < n; col+=8) { // doing 8 columns at a time
      REAL valA = *ptrToA;
      ptrToA++; // advance A to the next element in the row block
      REAL c0 = valA * ptrToB[col];
      REAL c1 = valA * ptrToB[col+1];
      REAL c2 = valA * ptrToB[col+2];
      REAL c3 = valA * ptrToB[col+3];
      REAL c4 = valA * ptrToB[col+4];
      REAL c5 = valA * ptrToB[col+5];
      REAL c6 = valA * ptrToB[col+6];
      REAL c7 = valA * ptrToB[col+7];

      // accumulate into c_i the product that needs to go into current row
      // block of C
      for(int j = 1; j < n; j++) {
        // each iter does a block of columns in the current row
        ptrToB += b_orig_n; // advance B to the next row
        valA = *ptrToA;
        ptrToA++; // advance A to the next element
        c0 += valA * ptrToB[col];
        c1 += valA * ptrToB[col+1];
        c2 += valA * ptrToB[col+2];
        c3 += valA * ptrToB[col+3];
        c4 += valA * ptrToB[col+4];
        c5 += valA * ptrToB[col+5];
        c6 += valA * ptrToB[col+6];
        c7 += valA * ptrToB[col+7];
      }
      // at this point ptrToB points to the bottom row
      // and ptrToA points to the end of a row
      //assert(ptrToA == &A[row*n + n]);
      //assert(ptrToB == &B[n*(n-1)]);
      ptrToB = B;  // put ptrToB back to the beginning of B
      ptrToA -= n; // move A back to beginning of the row
      *(C) = c0;
      *(C+1) = c1;
      *(C+2) = c2;
      *(C+3) = c3;
      *(C+4) = c4;
      *(C+5) = c5;
      *(C+6) = c6;
      *(C+7) = c7;
      //assert(C == &oldC[row*n + col]);
      C += 8;
      dummy += 8;
    }
    //assert(C == &oldC[row*n + n]);
    C += (c_orig_n - n);
    dummy += (c_orig_n - n);
    ptrToA += a_orig_n; // advance A to the next row
  }
  return;
}

// recursive parallel solution to matrix multiplication
void mm_dac_z(REAL *C, REAL *A, REAL *B, int n, bool add, int c_orig_n, int a_orig_n, int b_orig_n) {

  if(n == DAC_ARITH_BASECASE) {
    if(add) mm_additive_base(C, A, B, n, c_orig_n, a_orig_n, b_orig_n);
    else    mm_base(C, A, B, n, c_orig_n, a_orig_n, b_orig_n);
    return;
  }

  // partition each matrix into 4 sub matrices
  // each sub-matrix points to the start of the z pattern
  REAL *A11, *A12, *A21, *A22;
  REAL *B11, *B12, *B21, *B22;
  REAL *C11, *C12, *C21, *C22;
  Z_PARTITION(A, A11, A12, A21, A22, n, a_orig_n);
  Z_PARTITION(B, B11, B12, B21, B22, n, b_orig_n);
  Z_PARTITION(C, C11, C12, C21, C22, n, c_orig_n);

  // recrusively call the sub-matrices for evaluation in parallel
  cilk_spawn mm_dac_z(C11, A11, B11, n >> 1, add, c_orig_n, a_orig_n, b_orig_n);
  cilk_spawn mm_dac_z(C12, A11, B12, n >> 1, add, c_orig_n, a_orig_n, b_orig_n);
  cilk_spawn mm_dac_z(C21, A21, B11, n >> 1, add, c_orig_n, a_orig_n, b_orig_n);
  mm_dac_z(C22, A21, B12, n >> 1, add, c_orig_n, a_orig_n, b_orig_n);
  cilk_sync;

  cilk_spawn mm_dac_z(C11, A12, B21, n >> 1, true, c_orig_n, a_orig_n, b_orig_n);
  cilk_spawn mm_dac_z(C12, A12, B22, n >> 1, true, c_orig_n, a_orig_n, b_orig_n);
  cilk_spawn mm_dac_z(C21, A22, B21, n >> 1, true, c_orig_n, a_orig_n, b_orig_n);
  mm_dac_z(C22, A22, B22, n >> 1, true, c_orig_n, a_orig_n, b_orig_n);
  cilk_sync;

  return;
}

static void final_add(REAL *C11, REAL *C12, REAL *C21, REAL *C22,
                      REAL *M2, REAL *M5, REAL *T1, int n, int c_orig_n, int temp_orig_n) {

  if(n == DAC_ARITH_BASECASE) {
    int c_offset = c_orig_n - n;
    int temp_offset = temp_orig_n - n;
    for(int i=0; i<n; i++) { // this is only ok because we have square matrix
      for(int j=0; j<n; j+=8) {
        UNROLL(M5, m5); UNROLL(M2, m2); UNROLL(T1, t1);
        UNROLL_ADD_1(C11, m2);
        UNROLL_ADD_3(C12, m5, t1, m2);
        UNROLL(C22, c22);
        UNROLL_VAL_ADD_2(c22, m2, t1);
        UNROLL_VAL_SUB_SELF(C21, c22);
        UNROLL_ADD_2_ASSIGN(C22, c22, m5);
        M5 += 8; M2 += 8; T1 += 8;
        C11 += 8; C12 += 8; C21 += 8; C22 += 8;
      }
        M5 += temp_offset; M2 += temp_offset; T1 += temp_offset;
        C11 += c_offset; C12 += c_offset; C21 += c_offset; C22 += c_offset;
    }

    return;
  }

  REAL *C11_0, *C11_1, *C11_2, *C11_3;
  REAL *C12_0, *C12_1, *C12_2, *C12_3;
  REAL *C21_0, *C21_1, *C21_2, *C21_3;
  REAL *C22_0, *C22_1, *C22_2, *C22_3;
  Z_PARTITION(C11, C11_0, C11_1, C11_2, C11_3, n, c_orig_n);
  Z_PARTITION(C12, C12_0, C12_1, C12_2, C12_3, n, c_orig_n);
  Z_PARTITION(C21, C21_0, C21_1, C21_2, C21_3, n, c_orig_n);
  Z_PARTITION(C22, C22_0, C22_1, C22_2, C22_3, n, c_orig_n);

  REAL *M2_0, *M2_1, *M2_2, *M2_3;
  REAL *M5_0, *M5_1, *M5_2, *M5_3;
  REAL *T1_0, *T1_1, *T1_2, *T1_3;
  Z_PARTITION(M2, M2_0, M2_1, M2_2, M2_3, n, temp_orig_n);
  Z_PARTITION(M5, M5_0, M5_1, M5_2, M5_3, n, temp_orig_n);
  Z_PARTITION(T1, T1_0, T1_1, T1_2, T1_3, n, temp_orig_n);

  cilk_spawn final_add(C11_0, C12_0, C21_0, C22_0, M2_0, M5_0, T1_0, n >> 1, c_orig_n, temp_orig_n);
  cilk_spawn final_add(C11_1, C12_1, C21_1, C22_1, M2_1, M5_1, T1_1, n >> 1, c_orig_n, temp_orig_n);
  cilk_spawn final_add(C11_2, C12_2, C21_2, C22_2, M2_2, M5_2, T1_2, n >> 1, c_orig_n, temp_orig_n);
             final_add(C11_3, C12_3, C21_3, C22_3, M2_3, M5_3, T1_3, n >> 1, c_orig_n, temp_orig_n);
  cilk_sync;

  return;
}

#define Z_PARTITION_AND_DECL(S, n, orig_n) \
  REAL *S##_0, *S##_1, *S##_2, *S##_3; \
  S##_0 = &S[0]; \
  S##_1 = &S[n >> 1]; \
  S##_2 = &S[(n * orig_n) >> 1]; \
  S##_3 = &S[((n * orig_n) + n) >> 1];\

static void setup_add(REAL *S1, REAL *S2, REAL *S3, REAL *S4,
                      REAL *S5, REAL *S6, REAL *S7, REAL *S8,
                      REAL *A11, REAL *A12, REAL *A21, REAL *A22,
                      REAL *B11, REAL *B12, REAL *B21, REAL *B22, int n, int temp_orig_n, int a_orig_n, int b_orig_n) {

  if(n == DAC_ARITH_BASECASE) {
    int a_offset = a_orig_n - n;
    int b_offset = b_orig_n - n;
    int temp_offset = temp_orig_n - n;
    for(int i=0; i<n; i++) { // this is only ok because we have square matrix
      for(int j=0; j<n; j+=8) {
        UNROLL(A11, a11); UNROLL(A12, a12); UNROLL(A21, a21); UNROLL(A22, a22);
        UNROLL(B11, b11); UNROLL(B12, b12); UNROLL(B21, b21); UNROLL(B22, b22);
        UNROLL_ADD_2_ASSIGN(S1, a21, a22);
        UNROLL(S1, s1);
        UNROLL_SUB_2_ASSIGN(S2, s1, a11);
        UNROLL(S2, s2);
        UNROLL_SUB_2_ASSIGN(S4, a12, s2);
        UNROLL_SUB_2_ASSIGN(S3, a11, a21);
        UNROLL_SUB_2_ASSIGN(S5, b12, b11);
        UNROLL(S5, s5);
        UNROLL_SUB_2_ASSIGN(S6, b22, s5);
        UNROLL_SUB_2_ASSIGN(S7, b22, b12);
        UNROLL(S6, s6);
        UNROLL_SUB_2_ASSIGN(S8, s6, b21);

        A11 += 8; A12 += 8; A21 += 8; A22 += 8;
        B11 += 8; B12 += 8; B21 += 8; B22 += 8;
        S1 += 8; S2 += 8; S3 += 8; S4 += 8;
        S5 += 8; S6 += 8; S7 += 8; S8 += 8;
      }
        A11 += a_offset; A12 += a_offset; A21 += a_offset; A22 += a_offset;
        B11 += b_offset; B12 += b_offset; B21 += b_offset; B22 += b_offset;
        S1 += temp_offset; S2 += temp_offset; S3 += temp_offset; S4 += temp_offset;
        S5 += temp_offset; S6 += temp_offset; S7 += temp_offset; S8 += temp_offset;
    }
    return;
  }

  //Temp matricies need the "new n" from the calling strassen fuction
  //they shrink with each new lower recurssion level
  //A and B need orig_n as their size does not change
  Z_PARTITION_AND_DECL(S1, n, temp_orig_n); Z_PARTITION_AND_DECL(S2, n, temp_orig_n);
  Z_PARTITION_AND_DECL(S3, n, temp_orig_n); Z_PARTITION_AND_DECL(S4, n, temp_orig_n);
  Z_PARTITION_AND_DECL(S5, n, temp_orig_n); Z_PARTITION_AND_DECL(S6, n, temp_orig_n);
  Z_PARTITION_AND_DECL(S7, n, temp_orig_n); Z_PARTITION_AND_DECL(S8, n, temp_orig_n);
  Z_PARTITION_AND_DECL(A11, n, a_orig_n); Z_PARTITION_AND_DECL(A12, n, a_orig_n);
  Z_PARTITION_AND_DECL(A21, n, a_orig_n); Z_PARTITION_AND_DECL(A22, n, a_orig_n);
  Z_PARTITION_AND_DECL(B11, n, b_orig_n); Z_PARTITION_AND_DECL(B12, n, b_orig_n);
  Z_PARTITION_AND_DECL(B21, n, b_orig_n); Z_PARTITION_AND_DECL(B22, n, b_orig_n);

  cilk_spawn setup_add(S1_0,S2_0,S3_0,S4_0,S5_0,S6_0,S7_0,S8_0,
                       A11_0,A12_0,A21_0,A22_0,B11_0,B12_0,B21_0,B22_0,n>>1, temp_orig_n, a_orig_n, b_orig_n);
  cilk_spawn setup_add(S1_1,S2_1,S3_1,S4_1,S5_1,S6_1,S7_1,S8_1,
                       A11_1,A12_1,A21_1,A22_1,B11_1,B12_1,B21_1,B22_1,n>>1, temp_orig_n, a_orig_n, b_orig_n);
  cilk_spawn setup_add(S1_2,S2_2,S3_2,S4_2,S5_2,S6_2,S7_2,S8_2,
                       A11_2,A12_2,A21_2,A22_2,B11_2,B12_2,B21_2,B22_2,n>>1, temp_orig_n, a_orig_n, b_orig_n);
             setup_add(S1_3,S2_3,S3_3,S4_3,S5_3,S6_3,S7_3,S8_3,
                       A11_3,A12_3,A21_3,A22_3,B11_3,B12_3,B21_3,B22_3,n>>1, temp_orig_n, a_orig_n, b_orig_n);
  cilk_sync;
}

#define NUM_TEMP_M 11   // 11 temp matrices needed
/**
 * Perform C = A x B using strassen (except the base case)
 * n: size of current C, A, and B
 * It's assume that all C, A, B are using blocked z row-major layout
 **/
void strassen_z(REAL *C, REAL *A, REAL *B, int n, int c_orig_n, int a_orig_n, int b_orig_n) {

  if(n <= MATMUL_THRESH) {
    mm_dac_z(C, A, B, n, false, c_orig_n, a_orig_n, b_orig_n);
    return;
  }

  int new_n = (n >> 1);
  int subm_size = new_n * new_n; // in number of elements
  subm_size += (CACHE_LINE_SIZE / sizeof(REAL)); // avoid false sharing

  REAL *tmp = (REAL *) malloc(subm_size * sizeof(REAL) * NUM_TEMP_M);
  REAL *old_tmp = tmp;
  REAL *S1 = tmp; tmp += subm_size;
  REAL *S2 = tmp; tmp += subm_size;
  REAL *S3 = tmp; tmp += subm_size;
  REAL *S4 = tmp; tmp += subm_size;
  REAL *S5 = tmp; tmp += subm_size;
  REAL *S6 = tmp; tmp += subm_size;
  REAL *S7 = tmp; tmp += subm_size;
  REAL *S8 = tmp; tmp += subm_size;
  REAL *M2 = tmp; tmp += subm_size;
  REAL *M5 = tmp; tmp += subm_size;
  REAL *T1 = tmp;

  REAL *A11, *A12, *A21, *A22;
  REAL *B11, *B12, *B21, *B22;
  REAL *C11, *C12, *C21, *C22;
  Z_PARTITION(A, A11, A12, A21, A22, n, a_orig_n);
  Z_PARTITION(B, B11, B12, B21, B22, n, b_orig_n);
  Z_PARTITION(C, C11, C12, C21, C22, n, c_orig_n);

  /* The setup_add replace these; untested
  cilk_spawn add_matrix(S1, A21, A22, new_n);
  cilk_spawn sub_matrix(S3, A11, A21, new_n);
  cilk_spawn sub_matrix(S5, B12, B11, new_n);
             sub_matrix(S7, B22, B12, new_n);
  cilk_sync;

  sub_matrix(S2, S1, A11, new_n);
  sub_matrix(S4, A12, S2, new_n);
  sub_matrix(S6, B22, S5, new_n);
  sub_matrix(S8, S6, B21, new_n);
  */

  setup_add(S1,S2,S3,S4,S5,S6,S7,S8,A11,A12,A21,A22,B11,B12,B21,B22, new_n, new_n, a_orig_n, b_orig_n);

  cilk_spawn strassen_z(M2, A11, B11, new_n, new_n, a_orig_n, b_orig_n);  // P1, store in M2
  cilk_spawn strassen_z(M5, S1, S5, new_n, new_n, new_n, new_n);    // P5, store in M5
  cilk_spawn strassen_z(T1, S2, S6, new_n, new_n, new_n, new_n);    // P6, store in T1
  cilk_spawn strassen_z(C22, S3, S7, new_n, c_orig_n, new_n, new_n);   // P7, store in C22
  cilk_spawn strassen_z(C11, A12, B21, new_n, c_orig_n, a_orig_n, b_orig_n); // P2, store in C11
  cilk_spawn strassen_z(C12, S4, B22, new_n, c_orig_n, new_n, b_orig_n);  // P3, store in C12
             strassen_z(C21, A22, S8, new_n, c_orig_n, a_orig_n, new_n);  // P4, store in C21
  cilk_sync;

  final_add(C11, C12, C21, C22, M2, M5, T1, new_n, c_orig_n, new_n);

  free(old_tmp);
}

/*
 * Compare two matrices.  Print an error message if they differ by
 * more than EPSILON.
 * return 0 if no error, and return non-zero if error.
 */
static int compare_matrix(REAL *A, REAL *B, int n, int orig_n) {

  if(n == DAC_ARITH_BASECASE) {
    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        // compute the relative error c
        REAL c = A[orig_n*i + j] - B[orig_n*i + j];
        if (c < 0.0) c = -c;

        c = c / A[orig_n*i + j];
        if (c > EPSILON) {
          return -1;
        }
      }
    }
    return 0;
  }

  REAL *A11, *A12, *A21, *A22;
  REAL *B11, *B12, *B21, *B22;
  Z_PARTITION(A, A11, A12, A21, A22, n, orig_n);
  Z_PARTITION(B, B11, B12, B21, B22, n, orig_n);

  int r1 = compare_matrix(A11, B11, n >> 1, orig_n);
  int r2 = compare_matrix(A12, B12, n >> 1, orig_n);
  int r3 = compare_matrix(A21, B21, n >> 1, orig_n);
  int r4 = compare_matrix(A22, B22, n >> 1, orig_n);

  return (r1+r2+r3+r4);
}

int usage(void) {
  fprintf(stderr,
      "\nUsage: strassen_z [-n #] [-c]\n\n"
      "Multiplies two randomly generated n x n matrices. To check for\n"
      "correctness use -c\n");
  return 1;
}

const char *specifiers[] = {"-n", "-c", "-h", 0};
int opt_types[] = {INTARG, BOOLARG, BOOLARG, 0};

int main(int argc, char *argv[]) {

  REAL *A, *B, *C;
  int verify, help, n;

  /* standard benchmark options*/
  n = 2048;
  verify = 0;

  get_options(argc, argv, specifiers, opt_types, &n, &verify, &help);
  if (help || argc == 1) return usage();

  if((n & (n - 1)) != 0 || (n % 16) != 0) {
    fprintf(stderr, "%d: matrix size must be a power of 2"
            " and a multiple of %d\n", n, 16);
    return 1;
  }
  if(POWER < 3) {
    fprintf(stderr, "Must have at least base case of 8x8.\n");
    return 1;
  }

  A = (REAL *) malloc(n * n * sizeof(REAL));
  B = (REAL *) malloc(n * n * sizeof(REAL));
  C = (REAL *) malloc(n * n * sizeof(REAL));

  init_matrix(n, A, n);
  init_matrix(n, B, n);
  init_matrix(n, C, n);

  top_n = n;
  clockmark_t begin, end;

  __cilkrts_init();

#if TIMING_COUNT
  uint64_t elapsed_times[TIMING_COUNT];

  for(int i=0; i < TIMING_COUNT; i++) {
    __cilkrts_reset_timing();
    begin = ktiming_getmark();
    strassen_z(C, A, B, n, n, n, n);
    end = ktiming_getmark();
    elapsed_times[i] = ktiming_diff_usec(&begin, &end);
  }
  print_runtime(elapsed_times, TIMING_COUNT);
#else
  __cilkrts_reset_timing();
  begin = ktiming_getmark();
  strassen_z(C, A, B, n, n, n, n);
  end = ktiming_getmark();
  double elapsed = ktiming_diff_sec(&begin, &end);
  printf("Elapsed time in second: %f\n", elapsed);
#endif
  __cilkrts_accum_timing();

  if (verify) {
    printf("Checking results ... \n");
    REAL *C2 = (REAL *) malloc(n * n * sizeof(REAL));
    matrixmul(C2, A, B, n, n);
    verify = compare_matrix(C, C2, n, n);
    free(C2);
  }

  if(verify) {
    printf("WRONG RESULT!\n");
  } else {
    printf("\nCilk Example: strassen\n");
    printf("Options: n = %d\n\n", n);
  }

  free(A);
  free(B);
  free(C);

  return 0;
}
