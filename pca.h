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

#endif  // PCA_H_
