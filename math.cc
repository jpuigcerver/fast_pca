#include "math.h"

extern "C" {
  void saxpy_(int*, const float*, const float*, int*, float*, int*);
  void daxpy_(int*, const double*, const double*, int*, double*, int*);
  void sger_(int*, int*, float*, const float*, int*, const float*, int*,
             float*, int*);
  void dger_(int*, int*, double*, const double*, int*, const double*, int*,
             double*, int*);
  void ssyev_(char* jobz, char* uplo, int* n, float* a, int* lda, float*w,
             float* work, int* lwork, int* info);
  void dsyev_(char* jobz, char* uplo, int* n, double* a, int* lda, double*w,
             double* work, int* lwork, int* info);
}

template <> void axpy<float>(
    int n, float a, const float* const x, float* y) {
  int inc = 1;
  saxpy_(&n, &a, x, &inc, y, &inc);
}

template <> void axpy<double>(
    int n, double a, const double* const x, double* y) {
  int inc = 1;
  daxpy_(&n, &a, x, &inc, y, &inc);
}

template <> void ger<float>(
    int m, int n, float alpha, const float* x, const float* y, float* a) {
  int inc = 1;
  sger_(&m, &n, &alpha, x, &inc, y, &inc, a, &n);
}

template <> void ger<double>(
    int m, int n, double alpha, const double* x, const double* y, double* a) {
  int inc = 1;
  dger_(&m, &n, &alpha, x, &inc, y, &inc, a, &n);
}

template <> int syev<float>(int n, float* a, float* w) {
  char opt[2] = {'V', 'U'};
  int info = 0, lwork = -1;
  float wkopt = 0;
  ssyev_(opt, opt + 1, &n, a, &n, w, &wkopt, &lwork, &info);
  if (info != 0) { return info; }
  lwork = (int)wkopt;
  float* work = new float [lwork];
  ssyev_(opt, opt + 1, &n, a, &n, w, work, &lwork, &info);
  delete [] work;
  return info;
}

template <> int syev<double>(int n, double* a, double* w) {
  char opt[2] = {'V', 'U'};
  int info = 0, lwork = -1;
  double wkopt = 0;
  dsyev_(opt, opt + 1, &n, a, &n, w, &wkopt, &lwork, &info);
  if (info != 0) { return info; }
  lwork = (int)wkopt;
  double* work = new double [lwork];
  dsyev_(opt, opt + 1, &n, a, &n, w, work, &lwork, &info);
  delete [] work;
  return info;
}
