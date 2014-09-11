#ifndef PCA_H_
#define PCA_H_

#include <cmath>

#include "math.h"

template <typename real_t>
int pca(int dims, real_t* cov, real_t* w, real_t* s) {
  // compute standard deviation of each dimension
  for (int d = 0; d < dims; ++d) {
    s[d] = sqrt(cov[d * dims + d]);
  }
  // compute sorted eigenvectors and eigenvalues
  return eig<real_t>(dims, cov, w);
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
  gemv<real_t>('T', odim, idim, 1, eigvec, odim, B.data(), 1, 0, data, 1);
}

#endif  // PCA_H_
