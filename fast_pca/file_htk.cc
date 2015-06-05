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

#include "fast_pca/file_htk.h"
#include "fast_pca/endian.h"

// virtual
bool MatrixFile_HTK::copy_header_from(const MatrixFile& other) {
  if (other.format() != format_) return false;
  rows_ = other.rows();
  cols_ = other.cols();
  const MatrixFile_HTK* other_htk = static_cast<const MatrixFile_HTK*>(&other);
  nSamples_ = other_htk->nSamples_;
  sampPeriod_ = other_htk->sampPeriod_;
  sampSize_ = other_htk->sampSize_;
  parmKind_ = other_htk->parmKind_;
  return true;
}

// virtual
bool MatrixFile_HTK::read_header() {
  CHECK(file_);
  if (fread(&nSamples_, 4, 1, file_) != 1) return false;
  if (fread(&sampPeriod_, 4, 1, file_) != 1) return false;
  if (fread(&sampSize_, 2, 1, file_) != 1) return false;
  if (fread(&parmKind_, 2, 1, file_) != 1) return false;
  nSamples_   = be32toh(nSamples_);
  sampPeriod_ = be32toh(sampPeriod_);
  sampSize_   = be16toh(sampSize_);
  parmKind_   = be16toh(parmKind_);
  CHECK_FMT(
      sampSize_ % 4 == 0,
      "Currently, HTK format assumes that all input elements are floats, "
      "but the read sampSize (%d) is incompatible with this", sampSize_);
  rows_ = nSamples_;
  cols_ = sampSize_ / 4;
  return true;
}

// virtual
void MatrixFile_HTK::write_header() const {
  CHECK(file_);
  const uint32_t be_nSamples   = htobe32(rows_);
  const uint32_t be_sampPeriod = htobe32(sampPeriod_);
  const uint16_t be_sampSize   = htobe16(cols_ * 4);
  const uint16_t be_parmKind   = htobe16(parmKind_);
  fwrite(&be_nSamples, 4, 1, file_);
  fwrite(&be_sampPeriod, 4, 1, file_);
  fwrite(&be_sampSize, 2, 1, file_);
  fwrite(&be_parmKind, 2, 1, file_);
}

// virtual
int MatrixFile_HTK::read_block(int n, float* m) const {
  CHECK(file_);
  n = fread(m, 4, n, file_);
  for (int i = 0; i < n; ++i) {
    m[i] = beftoh(m[i]);
  }
  return n;
}

// virtual
int MatrixFile_HTK::read_block(int n, double* m) const {
  CHECK(file_);
  int i = 0;
  float f = 0;
  for (; i < n && fread(&f, 4, 1, file_) == 1; ++i) {
    m[i] = beftoh(f);
  }
  return i;
}

// virtual
void MatrixFile_HTK::write_block(int n, const float* m) const {
  CHECK(file_);
  for (int i = 0; i < n; ++i) {
    const float f = htobef(m[i]);
    fwrite(&f, 4, 1, file_);
  }
}

// virtual
void MatrixFile_HTK::write_block(int n, const double* m) const {
  CHECK(file_);
  for (int i = 0; i < n; ++i) {
    const float f = htobef(m[i]);
    fwrite(&f, 4, 1, file_);
  }
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_HTK>() {
  return new MatrixFile_HTK;
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_HTK>(FILE* file) {
  return new MatrixFile_HTK(file);
}
