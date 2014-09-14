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
#include <vector>

#include "fast_pca/math.h"

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

template <typename real_t>
int project(
    int n, int idim, int odim, const real_t* eigvec, const real_t* mean,
    const real_t* stdv, real_t* data) {
  if (idim < odim) { return -1; }
  // center input data (B = X - mean)
  vector<real_t> B(data, data + n * idim);
  // this can run in parallel
  for (int i = 0; i < n; ++i) {
    axpy<real_t>(idim, -1, mean, B.data() + i * idim);
  }
  // if the variances are given, input data is normalized
  if (stdv) {
    for (int i = 0; i < n * idim; ++i) {
      B[i] /= stdv[i % idim];
    }
  }
  gemm<real_t>(
      'N', 'N', n, odim, idim, 1, B.data(), n, eigvec, idim, 1, data, n);
  return 0;
}

template <typename real_t>
int project_single(
    int idim, int odim, const real_t* eigvec, const real_t* mean,
    const real_t* stdv, real_t* data) {
  if (idim < odim) { return -1; }
  // center input data (B = X - mean)
  vector<real_t> B(data, data + idim);
  axpy<real_t>(idim, -1, mean, B.data());
  // if the variances are given, input data is normalized
  if (stdv) {
    for (int i = 0; i < idim; ++i) {
      B[i] /= stdv[i];
    }
  }
  /*gemm<real_t>(
    'N', 'T', 1, odim, idim, 1, B.data(), idim, eigvec, odim, 1, data, odim);*/
  for (int od = 0; od < odim; ++od) {
    data[od] = 0;
    for (int id = 0; id < idim; ++id) {
      data[od] += B.data()[id] * eigvec[od * idim + id];
    }
  }
  return 0;
}

#endif  // FAST_PCA_PCA_H_
