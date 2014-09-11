#ifndef PARTIAL_H_
#define PARTIAL_H_

#include "math.h"
#include "file.h"

#include <cassert>
#include <numeric>
#include <vector>

using std::accumulate;
using std::vector;

// Compute sufficient variables to compute the covariance matrix
template <bool ascii, typename real_t>
int compute_partial(FILE* file, int dims, real_t* C, real_t* S) {
  assert(C); assert(S);
  vector<real_t> x(dims);
  int d = 0, n = 0;
  for (; (d = read_row<ascii, real_t>(file, dims, x.data())) == dims; ++n) {
    axpy<real_t>(dims, 1.0, x.data(), S);
    ger<real_t>(dims, dims, 1.0, x.data(), x.data(), C);
  }

  return (d == 0 || d == dims) ? n : -1;
}

template <typename real_t>
int reduce_partial(int mappers, int dims, int* n, real_t* C, real_t* S) {
  assert(n); assert(C); assert(S);
  const int N = accumulate(n, n + mappers, 0);
  for (int m = 1; m < mappers; ++m) {
    axpy<real_t>(dims, 1.0, S + m * dims, S);
    axpy<real_t>(dims * dims, 1.0, C + m * dims * dims, C);
  }
  return N;
}

#endif  // PARTIAL_H_
