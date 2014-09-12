#ifndef MATH_H_
#define MATH_H_

#include <algorithm>
#include <vector>

using std::sort;
using std::swap;
using std::vector;

// y += alpha * x
template <typename real_t>
void axpy(int n, real_t alpha, const real_t* x, real_t* y);

// A += alpha * x * y^T
template <typename real_t>
void ger(
    int m, int n, real_t alpha, const real_t* x, const real_t* y, real_t* A);

// compute eigenvectors and eigenvalues of a symmetric matrix
template <typename real_t>
int syev(int n, real_t* a, real_t* w);

// C = alpha * A * B + beta * C
template <typename real_t>
void gemm(
    char, char, int m, int n, int k, real_t alpha, const real_t* A, int lda,
    const real_t* B, int ldb, real_t beta, real_t* C, int ldc);

template <typename real_t>
void gemv();

// Compute eigenvalues and eigenvectors of the matrix m
// Important: Eigenvectors are stored as row vectors in the original matrix m
template <typename real_t>
int eig(int dims, real_t* m, real_t* w) {
  // Compute eigenvalues and eigenvectors
  const int info = syev<real_t>(dims, m, w);
  if (info != 0) { return info; }
  // Compute ordering of the eigenvalues
  vector<int> order(dims);
  for (int d = 0; d < dims; ++d) { order[d] = d; }
  sort(order.begin(), order.end(), [&w](int a, int b) -> bool {
      return w[a] > w[b];
    });
  // Reorder eigenvalues and eigenvectors
  for (int r = 0; r < dims / 2; ++r) {
    swap(w[r], w[order[r]]);
    for (int d = 0; d < dims; ++d) {
      real_t* x = m + r * dims + d;
      real_t* y = m + order[r] * dims + d;
      swap(*x, *y);
    }
  }
  return 0;
}

// axpy specializations for floats and doubles
template <> void axpy<float>(int, float, const float*, float*);
template <> void axpy<double>(int, double, const double*, double*);

// ger (outer product) specialization for float and doubles
template <> void ger<float>(
    int, int, float, const float*, const float*, float*);
template <> void ger<double>(
    int, int, double, const double*, const double*, double*);

// syev specializations for float and doubles
template <> int syev<float>(int n, float* a, float* w);
template <> int syev<double>(int n, double* a, double* w);

// gemm specializations for float and doubles
template <> void gemm<float>(
    char, char, int, int, int, float, const float*, int, const float*, int,
    float, float*, int);
template <> void gemm<double>(
    char, char, int, int, int, double, const double*, int, const double*, int,
    double, double*, int);

#endif  // MATH_H_
