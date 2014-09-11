#include "math.h"

extern "C" {
  void saxpy_(int*, const float*, const float*, int*, float*, int*);
  void daxpy_(int*, const double*, const double*, int*, double*, int*);
  void sger_(int*, int*, float*, const float*, int*, const float*, int*,
             float*, int*);
  void dger_(int*, int*, double*, const double*, int*, const double*, int*,
             double*, int*);
  void ssyev(char* jobz, char* uplo, int* n, float* a, int* lda, float*w,
             float* work, int* lwork, int* info);
  void dsyev(char* jobz, char* uplo, int* n, double* a, int* lda, double*w,
             double* work, int* lwork, int* info);
}

template <> void axpy<float>(
    int n, float a, const float* const x, int incx, float* y, int incy) {
  saxpy_(&n, &a, x, &incx, y, &incy);
}

template <> void axpy<double>(
    int n, double a, const double* const x, int incx, double* y, int incy) {
  daxpy_(&n, &a, x, &incx, y, &incy);
}

template <> void ger<float>(
    int m, int n, float alpha, const float* x, int incx, const float* y,
    int incy, float* a, int lda) {
  sger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}

template <> void ger<double>(
    int m, int n, double alpha, const double* x, int incx, const double* y,
    int incy, double* a, int lda) {
  dger_(&m, &n, &alpha, x, &incx, y, &incy, a, &lda);
}
