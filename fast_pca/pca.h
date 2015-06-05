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
// l -> (input)  leading dimension of matrix m
// m -> (input)  squared & symmetric matrix, (output) eigenvectors
// w -> (output) eigenvalues
template <typename real_t>
int eig(int n, int l, real_t* m, real_t* w) {
  // Compute eigenvalues and eigenvectors
  const int info = syev<real_t>(n, l, m, w);
  if (info != 0) { return info; }
  // Reorder eigenvalues and eigenvectors (ascending -> descending order)
  for (int r = 0; r < n / 2; ++r) {
    swap(w[r], w[n - r - 1]);
    for (int d = 0; d < n; ++d) {
      swap(m[r * l + d], m[(n - r - 1) * l + d]);
    }
  }
  return 0;
}

// n -> (input)  number of data samples
// p -> (input)  input data dimension
// q -> (input)  output data dimension
// r -> (input)  exclude these number of first/last dimensions from projection
// m -> (input)  mean of the input data for each dimension
// s -> (input)  standard deviation of the input data for each dimension
// v -> (input)  eigenvectors of the zero-mean covariance of the input data
// x -> (input)  original data,
//      (output) mean-centered and (optionally) normalized data
// z -> (output) projected data
template <typename real_t>
int project(
    int n, int p, int q, int r, const real_t* v, const real_t* m,
    const real_t* s, real_t* x, real_t* z) {
  if (p < q || !x || !z) { return -1; }
  // convert input data to zero-mean
  // TODO(jpuigcerver): this can run in parallel
  for (int i = 0; i < n; ++i) { axpy<real_t>(p, -1, m, x + i * p); }
  // if the variances are given, input data is normalized
  if (s) {
    // TODO(jpuigcerver): this can run in parallel
    for (int i = 0; i < n * p; ++i) {
      // just for safety, if the variance is small,
      // do not normalize data in that dimension
      if (s[i % p] > 1E-6) {
        x[i] /= s[i % p];
      }
    }
  }
  // Copy non-projected dimensions
  // TODO(jpuigcerver): this can run in parallel
  for (int i = 0; r != 0 && i < n; ++i) {
    const real_t* xi = r > 0 ? x + i * p : x + (i + 1) * p + r;
    real_t* zi = r > 0 ? z + i * q : z + (i + 1) * q + r;
    memcpy(zi, xi, sizeof(real_t) * abs(r));
  }
  // Effective sizes (p: input, q: output) of the projected data
  const int eff_p = p - abs(r);
  const int eff_q = q - abs(r);
  if (eff_p > 0 && eff_q > 0) {
    // Offset to the input/output data, so the excluded dimensions are not
    // projected
    real_t* x_offset = r <= 0 ? x : x + r;
    real_t* z_offset = r <= 0 ? z : z + r;
    gemm<real_t>(
        'N', 'T', n, eff_q, eff_p, 1, x_offset, p, v, eff_p, 0, z_offset, q);
  }
  return 0;
}

#endif  // FAST_PCA_PCA_H_
