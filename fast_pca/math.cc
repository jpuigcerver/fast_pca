/*
  The MIT License (MIT)

  Copyright (c) 2014 Joan Puigcerver

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include "fast_pca/math.h"

extern "C" {
  void saxpy_(int*, const float*, const float*, int*, float*, int*);
  void daxpy_(int*, const double*, const double*, int*, double*, int*);
  void sger_(int*, int*, float*, const float*, int*, const float*, int*,
             float*, int*);
  void dger_(int*, int*, double*, const double*, int*, const double*, int*,
             double*, int*);
  void ssyev_(char*, char*, int*, float*, int*, float*, float*, int*, int*);
  void dsyev_(char*, char*, int*, double*, int*, double*, double*, int*, int*);
  void sgemm_(char*, char*, int*, int*, int*, float*, const float*,
              int*, const float*, int*, float*, float*, int*);
  void dgemm_(char*, char*, int*, int*, int*, double*, const double*,
              int*, const double*, int*, double*, double*, int*);
}

template <> void axpy<float>(
    int n, float a, const float* const x, float* y) {
  int inc = 1;
  saxpy_(&n, &a, x, &inc, y, &inc);
}

template <> void axpy<double>(
    int n, double a, const double* const x, double* y) {
  int inc = 1;
  daxpy_(&n, &a, x, &inc, y, &inc);
}

template <> void ger<float>(
    int m, int n, float alpha, const float* x, const float* y, float* a) {
  int inc = 1;
  sger_(&m, &n, &alpha, x, &inc, y, &inc, a, &n);
}

template <> void ger<double>(
    int m, int n, double alpha, const double* x, const double* y, double* a) {
  int inc = 1;
  dger_(&m, &n, &alpha, x, &inc, y, &inc, a, &n);
}

template <> int syev<float>(int n, float* a, float* w) {
  char opt[2] = {'V', 'U'};
  int info = 0, lwork = -1;
  float wkopt = 0;
  ssyev_(opt, opt + 1, &n, a, &n, w, &wkopt, &lwork, &info);
  if (info != 0) { return info; }
  lwork = static_cast<int>(wkopt);
  float* work = new float[lwork];
  ssyev_(opt, opt + 1, &n, a, &n, w, work, &lwork, &info);
  delete [] work;
  return info;
}

template <> int syev<double>(int n, double* a, double* w) {
  char opt[2] = {'V', 'U'};
  int info = 0, lwork = -1;
  double wkopt = 0;
  dsyev_(opt, opt + 1, &n, a, &n, w, &wkopt, &lwork, &info);
  if (info != 0) { return info; }
  lwork = static_cast<int>(wkopt);
  double* work = new double[lwork];
  dsyev_(opt, opt + 1, &n, a, &n, w, work, &lwork, &info);
  delete [] work;
  return info;
}

void gemm_op(char* opA, char* opB) {
  char TA, TB;
  // determine op(B) in col-major order
  if (*opA == 'T') { TB = 'T'; }
  else if (*opA == 'C') { TB = 'C'; }
  else { TB = 'N'; }
  // determine op(A) in col-major order
  if (*opB == 'T') { TA = 'T'; }
  else if (*opB == 'C') { TA = 'C'; }
  else { TA = 'N'; }
  *opA = TA;
  *opB = TB;
}

template <> void gemm<float>(
    char opA, char opB, int m, int n, int k, float alpha, const float* A,
    int lda, const float* B, int ldb, float beta, float* C, int ldc) {
  gemm_op(&opA, &opB);
  sgemm_(&opA, &opB, &n, &m, &k, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
}

template <> void gemm<double>(
    char opA, char opB, int m, int n, int k, double alpha, const double* A,
    int lda, const double* B, int ldb, double beta, double* C, int ldc) {
  gemm_op(&opA, &opB);
  dgemm_(&opA, &opB, &n, &m, &k, &alpha, B, &ldb, A, &lda, &beta, C, &ldc);
}
