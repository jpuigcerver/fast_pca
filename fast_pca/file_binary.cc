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

#include "fast_pca/file_binary.h"

#include <cstdio>

// virtual
int MatrixFile_Binary::read_block(int n, float* m) const {
  CHECK(file_);
  return fread(m, sizeof(float), n, file_);
}

// virtual
int MatrixFile_Binary::read_block(int n, double* m) const {
  CHECK(file_);
  return fread(m, sizeof(double), n, file_);
}

// virtual
void MatrixFile_Binary::write_block(int n, const float* m) const {
  CHECK(file_);
  fwrite(m, sizeof(float), n, file_);
}

// virtual
void MatrixFile_Binary::write_block(int n, const double* m) const {
  CHECK(file_);
  fwrite(m, sizeof(double), n, file_);
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_BINARY>() {
  return new MatrixFile_Binary;
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_BINARY>(FILE* file) {
  return new MatrixFile_Binary(file);
}
