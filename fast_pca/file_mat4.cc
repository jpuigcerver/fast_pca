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

#include "fast_pca/file_mat4.h"
#include "fast_pca/endian.h"

// read from a file, using a given type TF, swap the bytes optionally, and
// then write to the buffer with type TT
template <typename TF, typename TT, bool swap>
int read_swap_cast_block(FILE* f, int n, TT* m) {
  TF t = 0;
  int i = 0;
  for (; i < n && fread(&t, sizeof(TF), 1, f) == 1; ++i) {
    if (swap) swap_bytes<sizeof(TF)>(&t);
    m[i] = t;
  }
  return i;
}

// write a block of data to a file, data is first casted to the output
// type and, optionally, bytes are swapped.
template <typename TF, typename TT, bool swap>
void write_swap_cast_block(FILE* f, int n, const TF* m) {
  for (int i = 0; i < n; ++i) {
    TT t = m[i];
    if (swap) swap_bytes<sizeof(TT)>(&t);
    fwrite(&t, sizeof(TT), 1, f);
  }
}

// virtual
bool MatrixFile_MAT4::copy_header_from(const MatrixFile& other) {
  if (other.format() != format_) return false;
  rows_ = other.rows();
  cols_ = other.cols();
  const MatrixFile_MAT4* other_mat4 =
      static_cast<const MatrixFile_MAT4*>(&other);
  name_ = other_mat4->name_;
  mopt_ = other_mat4->mopt_;
  prec_ = other_mat4->prec_;
  swap_ = other_mat4->swap_;
  return true;
}

// virtual
bool MatrixFile_MAT4::read_header() {
  CHECK(file_);
  uint32_t mrows;
  uint32_t ncols;
  uint32_t imagf;
  uint32_t namlen;
  if (fread(&mopt_, 4, 1, file_) != 1) return false;
  if (fread(&mrows, 4, 1, file_) != 1) return false;
  if (fread(&ncols, 4, 1, file_) != 1) return false;
  if (fread(&imagf, 4, 1, file_) != 1) return false;
  if (fread(&namlen, 4, 1, file_) != 1) return false;
  swap_ = false;
  if ((mopt_ == 0 && is_big_endian()) || mopt_ > 9999) {
    swap_bytes<4>(&mopt_);
    swap_bytes<4>(&mrows);
    swap_bytes<4>(&ncols);
    swap_bytes<4>(&imagf);
    swap_bytes<4>(&namlen);
    swap_ = true;
  }
  CHECK_FMT(
      mopt_ <= 9999, "Unsupported `mopt' field in MAT-v4 header (%u)!", mopt_);
  CHECK_FMT(
      imagf == 0, "Unsupported `imgf' field in MAT-v4 header (%u)!", imagf);
  if (namlen > 0) {
    char* buff = new char[namlen];
    if (fread(buff, 1, namlen, file_) != namlen) {
      delete [] buff;
      return false;
    }
    buff[namlen-1] = '\0';
    name_.assign(buff, namlen - 1);
  }
  CHECK_MSG(
      mopt_ % 10 == 0, "Only full matrix are supported by this "
      "implementation of the MAT-v4 format");
  // Check architecture that generated the file
  const char mach  = (mopt_ / 1000) % 10;
  CHECK_MSG(
      mach == 0 || mach == 1,
      "Only IEEE-754 Little Endian or IEEE-754 Big Endian matrix are "
      "supported by this implementation of the MAT-v4 format");
  // Col-major vs. row-major ordering
  order_ = (mopt_ / 100) % 10;
  if (order_ == 0) {
    // If MAT-V4 files are stored in col-major order, we must transpose the
    // matrix, since fast_pca assumes that the data samples are stored in
    // each row...
    rows_ = ncols;
    cols_ = mrows;
  } else {
    rows_ = mrows;
    cols_ = ncols;
  }
  // base data type: double, float, int, etc.
  prec_ = (mopt_ / 10) % 10;
  return true;
}

// virtual
void MatrixFile_MAT4::write_header() const {
  CHECK(file_);
  const uint32_t mopt = (is_big_endian() ? 1000 : 0) + \
      order_ * 100 + prec_ * 10;
  const uint32_t imagf = 0;
  const uint32_t namlen = name_.length() + 1;
  const uint32_t mrows = order_ == 0 ? cols_ : rows_;
  const uint32_t ncols = order_ == 0 ? rows_ : cols_;
  fwrite(&mopt,  4, 1, file_);
  fwrite(&mrows, 4, 1, file_);
  fwrite(&ncols, 4, 1, file_);
  fwrite(&imagf, 4, 1, file_);
  fwrite(&namlen, 4, 1, file_);
  fwrite(name_.c_str(), 1, namlen, file_);
}

template <typename T>
int MatrixFile_MAT4::read_block(int n, T* m) const {
  CHECK(file_);
  if (prec_ == type2prec<T>::prec) {
    n = fread(m, sizeof(T), n, file_);
    if (swap_) {
      for (int i = 0; i < n; ++i) swap_bytes<sizeof(T)>(m + i);
    }
    return n;
  } else if (prec_ == 0) {
    if (swap_)
      return read_swap_cast_block<double, T, true>(file_, n, m);
    else
      return read_swap_cast_block<double, T, false>(file_, n, m);
  } else if (prec_ == 1) {
    if (swap_)
      return read_swap_cast_block<float, T, true>(file_, n, m);
    else
      return read_swap_cast_block<float, T, false>(file_, n, m);
  } else {
    // TODO(jpuigcerver) Support additional casting
    ERROR_FMT(
        "With MAT-v4 cannot read from type %d to %d", prec_,
        type2prec<T>::prec);
  }
}

template <typename T>
void MatrixFile_MAT4::write_block(int n, const T* m) const {
  CHECK(file_);
  if (prec_ == type2prec<T>::prec) {
    if (swap_)
      write_swap_cast_block<T, T, true>(file_, n, m);
    else
      fwrite(m, sizeof(T), n, file_);
  } else if (prec_ == 0) {
    if (swap_)
      write_swap_cast_block<T, double, true>(file_, n, m);
    else
      write_swap_cast_block<T, double, false>(file_, n, m);
  } else if (prec_ == 1) {
    if (swap_)
      write_swap_cast_block<T, float, true>(file_, n, m);
    else
      write_swap_cast_block<T, float, false>(file_, n, m);
  } else {
    // TODO(jpuigcerver) Support additional casting
    ERROR_FMT(
        "With MAT-v4 cannot write from type %d to %d", type2prec<T>::prec,
        prec_);
  }
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_MAT4>() {
  return new MatrixFile_MAT4;
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_MAT4>(FILE* file) {
  return new MatrixFile_MAT4(file);
}

// Some static initialization
template <>
uint8_t MatrixFile_MAT4::type2prec<double>::prec   = 0;
template <>
uint8_t MatrixFile_MAT4::type2prec<float>::prec    = 1;
template <>
uint8_t MatrixFile_MAT4::type2prec<int32_t>::prec  = 2;
template <>
uint8_t MatrixFile_MAT4::type2prec<int16_t>::prec  = 3;
template <>
uint8_t MatrixFile_MAT4::type2prec<uint16_t>::prec = 4;
template <>
uint8_t MatrixFile_MAT4::type2prec<uint8_t>::prec  = 5;

// Instantiate read_block for different data types
template
int MatrixFile_MAT4::read_block(int, float*) const;
template
int MatrixFile_MAT4::read_block(int, double*) const;
template
int MatrixFile_MAT4::read_block(int, int32_t*) const;
template
int MatrixFile_MAT4::read_block(int, int16_t*) const;
template
int MatrixFile_MAT4::read_block(int, uint16_t*) const;
template
int MatrixFile_MAT4::read_block(int, uint8_t*) const;

// Instantiate write_block for different data types
template
void MatrixFile_MAT4::write_block(int, const float*) const;
template
void MatrixFile_MAT4::write_block(int, const double*) const;
template
void MatrixFile_MAT4::write_block(int, const int32_t*) const;
template
void MatrixFile_MAT4::write_block(int, const int16_t*) const;
template
void MatrixFile_MAT4::write_block(int, const uint16_t*) const;
template
void MatrixFile_MAT4::write_block(int, const uint8_t*) const;
