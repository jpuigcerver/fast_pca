#ifndef MATH_H_
#define MATH_H_

template <typename real_t>
void axpy(int n, real_t a, const real_t* x, int incx, real_t* y, int incy);

template <typename real_t>
void ger(int m, int n, real_t a, const real_t* x, int incx, const real_t* y,
         int incy, real_t* A, int lda);

template <typename real_t>
void syev();


template <> void axpy<float>(int, float, const float*, int, float*, int);
template <> void axpy<double>(int, double, const double*, int, double*, int);

template <> void ger<float>(
    int, int, float, const float*, int, const float*, int, float*, int);
template <> void ger<double>(
    int, int, double, const double*, int, const double*, int, double*, int);

#endif  // MATH_H_
