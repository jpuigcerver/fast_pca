/*
  The MIT License (MIT)

  Copyright (c) 2014,2015 Joan Puigcerver

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

#ifndef FAST_PCA_MATH_H_
#define FAST_PCA_MATH_H_

// y += alpha * x
template <typename real_t>
void axpy(int n, real_t alpha, const real_t* x, real_t* y);

// A += alpha * x * y^T
template <typename real_t>
void ger(
    int m, int n, real_t alpha, const real_t* x, const real_t* y, real_t* A);

// compute eigenvectors and eigenvalues of a symmetric matrix
template <typename real_t>
int syev(int n, int lda, real_t* a, real_t* w);

// C = alpha * A * B + beta * C
template <typename real_t>
void gemm(
    char, char, int m, int n, int k, real_t alpha, const real_t* A, int lda,
    const real_t* B, int ldb, real_t beta, real_t* C, int ldc);

// y = alpha * A * x + beta * y
template <typename real_t>
void gemv(
    char, int m, int n, real_t alpha, const real_t* A, int lda,
    const real_t* x, int incx, real_t beta, real_t* y, int incy);

// axpy specializations for floats and doubles
template <> void axpy<float>(int, float, const float*, float*);
template <> void axpy<double>(int, double, const double*, double*);

// ger (outer product) specialization for float and doubles
template <> void ger<float>(
    int, int, float, const float*, const float*, float*);
template <> void ger<double>(
    int, int, double, const double*, const double*, double*);

// syev specializations for float and doubles
template <> int syev<float>(int, int, float*, float*);
template <> int syev<double>(int, int, double*, double*);

// gemm specializations for float and doubles
template <> void gemm<float>(
    char, char, int, int, int, float, const float*, int, const float*, int,
    float, float*, int);
template <> void gemm<double>(
    char, char, int, int, int, double, const double*, int, const double*, int,
    double, double*, int);

// gemv specializations for float and doubles
template <> void gemv<float>(
    char, int, int, float, const float*, int, const float*, int, float,
    float*, int);
template <> void gemv<double>(
    char, int, int, double, const double*, int, const double*, int, double,
    double*, int);

#endif  // FAST_PCA_MATH_H_
