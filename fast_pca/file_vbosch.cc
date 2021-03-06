/*
  The MIT License (MIT)

  Copyright (c) 2014,2015 Joan Puigcerver

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

#include "fast_pca/file_vbosch.h"

// virtual
bool MatrixFile_VBosch::copy_header_from(const MatrixFile& other) {
  if (other.format() != format_) return false;
  rows_ = other.rows();
  cols_ = other.cols();
  return true;
}

// virtual
bool MatrixFile_VBosch::read_header() {
  CHECK(file_);
  return (fscanf(file_, "%d %d", &rows_, &cols_) == 2 &&
          rows_ >= 0 && cols_ >= 0);
}

// virtual
void MatrixFile_VBosch::write_header() const {
  CHECK(file_);
  fprintf(file_, "%d %d\n", rows_, cols_);
}

// virtual
int MatrixFile_VBosch::read_block(int n, float* m) const {
  CHECK(file_);
  int i = 0;
  for (; i < n && fscanf(file_, "%f", m + i) == 1; ++i) { }
  return i;
}

// virtual
int MatrixFile_VBosch::read_block(int n, double* m) const {
  CHECK(file_);
  int i = 0;
  for (; i < n && fscanf(file_, "%lf", m + i) == 1; ++i) { }
  return i;
}

// virtual
void MatrixFile_VBosch::write_block(int n, const float* m) const {
  CHECK(file_);
  for (int i = 0; i < n; ++i) {
    if ((i + 1) % cols_ == 0)
      fprintf(file_, "%.16g\n", m[i]);
    else
      fprintf(file_, "%.16g ", m[i]);
  }
}

// virtual
void MatrixFile_VBosch::write_block(int n, const double* m) const {
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
MatrixFile* MatrixFile::Create<FMT_VBOSCH>() {
  return new MatrixFile_VBosch;
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_VBOSCH>(FILE* file) {
  return new MatrixFile_VBosch(file);
}
