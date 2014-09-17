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

#include <cstdio>

#include "fast_pca/file.h"

// Plain ASCII format
template <> int read_matrix_header<FMT_ASCII>(
    FILE* file, string* name, int* rows, int* cols) {
  return 0;
}

template <> void write_matrix_header<FMT_ASCII>(
    FILE* file, const string& name, int rows, int cols) {
}

template <> int read_block<FMT_ASCII, float>(FILE* file, int n, float* m) {
  int i = 0;
  for (; i < n && fscanf(file, "%f", m + i) == 1; ++i) { }
  return i;
}

template <> int read_block<FMT_ASCII, double>(FILE* file, int n, double* m) {
  int i = 0;
  for (; i < n && fscanf(file, "%lf", m + i) == 1; ++i) { }
  return i;
}

template <> void write_block<FMT_ASCII, float>(
    FILE* file, int n, const float* m) {
  for (int i = 0; i < n - 1; ++i) fprintf(file, "%.12g ", m[i]);
  fprintf(file, "%.12g\n", m[n - 1]);
}

template <> void write_block<FMT_ASCII, double>(
    FILE* file, int n, const double* m) {
  for (int i = 0; i < n - 1; ++i) fprintf(file, "%.12g ", m[i]);
  fprintf(file, "%.12g\n", m[n - 1]);
}

// Plain Binary format
template <> int read_matrix_header<FMT_BINARY>(
    FILE* file, string* name, int* rows, int* cols) {
  return 0;
}

template <> void write_matrix_header<FMT_BINARY>(
    FILE* file, const string& name, int rows, int cols) {
}

template <> int read_block<FMT_BINARY, float>(FILE* file, int n, float* m) {
  return fread(m, sizeof(float), n, file);
}

template <> int read_block<FMT_BINARY, double>(FILE* file, int n, double* m) {
  return fread(m, sizeof(double), n, file);
}

template <> void write_block<FMT_BINARY, float>(
    FILE* file, int n, const float* m) {
  fwrite(m, sizeof(float), n, file);
}

template <> void write_block<FMT_BINARY, double>(
    FILE* file, int n, const double* m) {
  fwrite(m, sizeof(double), n, file);
}
