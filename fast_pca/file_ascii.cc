/*
  The MIT License (MIT)

  Copyright (c) 2015 Joan Puigcerver

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

#include "fast_pca/file_ascii.h"

#include <cstdio>

int MatrixFile_ASCII::read_block(int n, float* m) const {
  CHECK(file_);
  int i = 0;
  for (; i < n && fscanf(file_, "%f", m + i) == 1; ++i) { }
  return i;
}

int MatrixFile_ASCII::read_block(int n, double* m) const {
  CHECK(file_);
  int i = 0;
  for (; i < n && fscanf(file_, "%lf", m + i) == 1; ++i) { }
  return i;
}

void MatrixFile_ASCII::write_block(int n, const float* m) const {
  CHECK(file_);
  for (int i = 0; i < n; ++i) {
    if ((i + 1) % cols_ == 0)
      fprintf(file_, "%.16g\n", m[i]);
    else
      fprintf(file_, "%.16g ", m[i]);
  }
}

void MatrixFile_ASCII::write_block(int n, const double* m) const {
  CHECK(file_);
  for (int i = 0; i < n; ++i) {
    if ((i + 1) % cols_ == 0)
      fprintf(file_, "%.16g\n", m[i]);
    else
      fprintf(file_, "%.16g ", m[i]);
  }
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_ASCII>() {
  return new MatrixFile_ASCII;
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_ASCII>(FILE* file) {
  return new MatrixFile_ASCII(file);
}
