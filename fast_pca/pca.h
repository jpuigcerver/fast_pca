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

#ifndef FAST_PCA_PCA_H_
#define FAST_PCA_PCA_H_

#include <cmath>
#include <cstring>
#include <vector>

#include "fast_pca/math.h"
#include "fast_pca/file.h"

// n -> (input) number of data samples
// d -> (input) data dimension
// m -> (input) sum of each dimension, (output) mean of each dimension
// c -> (input) X' * X, (output) eigenvectors
// s -> (output) standard deviation of each dimension
// w -> (output) eigenvalues
template <typename real_t>
int pca(int n, int d, real_t* m, real_t* c, real_t* s, real_t* w) {
  // compute means
  for (int i = 0; i < d; ++i) { m[i] /= n; }
  // compute covariance matrix
  for (int i = 0; i < d; ++i) {
    for (int j = 0; j < d; ++j) {
      c[i * d + j] = (c[i * d + j] - n * m[i] * m[j]) / (n - 1);
    }
  }
  // compute standard deviation of each dimension
  for (int i = 0; i < d; ++i) {
    s[i] = sqrt(c[i * d + i]);
  }
  // compute sorted eigenvectors and eigenvalues
  return eig<real_t>(d, c, w);
}

// n -> (input) number of data samples
// p -> (input) input data dimension
// q -> (input) output data dimension
// e -> (input) eigenvectors of the zero-mean covariance of the input data
// m -> (input) mean of the input data for each dimension
// s -> (input) standard deviation of the input data for each dimension
// x -> (input/output) data
template <typename real_t>
int project(
    int n, int p, int q, const real_t* e, const real_t* m, const real_t* s,
    real_t* x) {
  if (p < q) { return -1; }
  // allocate space for zero-mean data
  vector<real_t> B(x, x + n * p);
  // convert input data to zero-mean
  // TODO(jpuigcerver): this can run in parallel
  for (int i = 0; i < n; ++i) { axpy<real_t>(p, -1, m, B.data() + i * p); }
  // if the variances are given, input data is normalized
  if (s) {
    // TODO(jpuigcerver): this can run in parallel
    for (int i = 0; i < n * p; ++i) { B[i] /= s[i % p]; }
  }
  gemm<real_t>('N', 'T', n, q, p, 1, B.data(), p, e, p, 0, x, q);
  return 0;
}

#endif  // FAST_PCA_PCA_H_
