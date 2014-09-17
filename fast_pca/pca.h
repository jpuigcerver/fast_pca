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

#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "fast_pca/math.h"

using std::sort;
using std::swap;
using std::vector;

// Compute eigenvalues and eigenvectors of the matrix m
// n -> (input)  number of dimensions
// m -> (input)  squared matrix, (output) eigenvectors
// w -> (output) eigenvalues
template <typename real_t>
int eig(int n, real_t* m, real_t* w) {
  // Compute eigenvalues and eigenvectors
  const int info = syev<real_t>(n, m, w);
  if (info != 0) { return info; }
  // Compute ordering of the eigenvalues
  vector<int> order(n);
  for (int d = 0; d < n; ++d) { order[d] = d; }
  sort(order.begin(), order.end(), [&w](int a, int b) -> bool {
      return w[a] > w[b];
    });
  // Reorder eigenvalues and eigenvectors
  for (int r = 0; r < n / 2; ++r) {
    swap(w[r], w[order[r]]);
    for (int d = 0; d < n; ++d) {
      real_t* x = m + r * n + d;
      real_t* y = m + order[r] * n + d;
      swap(*x, *y);
    }
  }
  return 0;
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
    for (int i = 0; i < n * p; ++i) {
      // just for safety, if the variance is small,
      // do not normalize data in that dimension
      if (s[i % p] > 1E-6) {
        B[i] /= s[i % p];
      }
    }
  }
  gemm<real_t>('N', 'T', n, q, p, 1, B.data(), p, e, p, 0, x, q);
  return 0;
}

#endif  // FAST_PCA_PCA_H_
