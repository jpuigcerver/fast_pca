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

#include "fast_pca/file.h"
#include "fast_pca/logging.h"

// Determine whether the machine is bigendian or not
// NOTE: The output of this function is known at compile time, a good
// compiler should know this (i.e. GCC does know it at compile time).
inline bool is_big_endian() {
  union {
    uint32_t i;
    char c[4];
  } bint = {0x01020304};
  return bint.c[0] == 1;
}

// Convert float number in big-endian to the host endianness
// NOTE: A clever compiler should optimize this heavily
inline float beftoh(float f) {
  if (is_big_endian()) return f;
  float f2 = 0.0;
  char* cf  = reinterpret_cast<char*>(&f);
  char* cf2 = reinterpret_cast<char*>(&f2);
  for (size_t i = 0; i < sizeof(float) / 2; ++i)
    cf[i] = cf2[sizeof(float) - i - 1];
  return f2;
}

// Convert float in host endianness to big-endian
// NOTE: A clever compiler should optimize this heavily
inline float htobef(float f) {
  if (is_big_endian()) return f;
  float f2 = 0.0;
  char* cf  = reinterpret_cast<char*>(&f);
  char* cf2 = reinterpret_cast<char*>(&f2);
  for (size_t i = 0; i < sizeof(float) / 2; ++i)
    cf[i] = cf2[sizeof(float) - i - 1];
  return f2;
}

template <> int read_matrix_header<FMT_PRHLT_HTK>(
    FILE* file, string* name, int* rows, int* cols) {
  uint32_t nSamples = 0, sampPeriod = 0;
  uint16_t sampSize = 0, parmKind = 0;
  if (fread(&nSamples, 4, 1, file) != 1) return 1;
  if (fread(&sampPeriod, 4, 1, file) != 1) return 1;
  if (fread(&sampSize, 2, 1, file) != 1) return 1;
  if (fread(&parmKind, 2, 1, file) != 1) return 1;
  *rows = static_cast<int>(be32toh(nSamples));
  *cols = static_cast<int>(be16toh(sampSize) / 4);
  parmKind = be16toh(parmKind);
  sampPeriod = be32toh(sampPeriod);
  CHECK_MSG(
      parmKind == 9, "HTK format only supports the USER parameter kind!");
  if (sampPeriod != 100000) {
    WARN_FMT("HTK format will save the data using a sampPeriod = 10Kus, "
         "but your data was originally sampled at %gKus...",
         sampPeriod * 0.1);
  }
  return 0;
}

template <> void write_matrix_header<FMT_PRHLT_HTK>(
    FILE* file, const string& name, int rows, int cols) {
  uint32_t nSamples = rows, sampPeriod = 1000000;
  uint16_t sampSize = cols * 4, parmKind = 9;
  nSamples   = htobe32(nSamples);
  sampPeriod = htobe32(sampPeriod);
  sampSize   = htobe16(sampSize);
  parmKind   = htobe16(parmKind);
  fwrite(&nSamples, 4, 1, file);
  fwrite(&sampPeriod, 4, 1, file);
  fwrite(&sampSize, 2, 1, file);
  fwrite(&parmKind, 2, 1, file);
}

template <>
int read_block<FMT_PRHLT_HTK, float>(FILE* file, int n, float* m) {
  int i = 0;
  for (; i < n && fread(m + i, 4, 1, file) == 1; ++i) {
    m[i] = beftoh(m[i]);
  }
  return i;
}

template <>
int read_block<FMT_PRHLT_HTK, double>(FILE* file, int n, double* m) {
  int i = 0;
  float f = 0;
  for (; i < n && fread(&f, 4, 1, file) == 1; ++i) {
    m[i] = beftoh(f);
  }
  return i;
}

template <>
void write_block<FMT_PRHLT_HTK, float>(FILE* file, int n, const float* m) {
  for (int i = 0; i < n; ++i) {
    const float f = htobef(m[i]);
    fwrite(&f, 4, 1, file);
  }
}

template <>
void write_block<FMT_PRHLT_HTK, double>(FILE* file, int n, const double* m) {
  for (int i = 0; i < n; ++i) {
    const double d = htobef(m[i]);
    fwrite(&d, 4, 1, file);
  }
}

template <>
void write_matrix<FMT_PRHLT_HTK, float>(
    FILE* file, int rows, int cols, const float* m) {
  for (int r = 0; r < rows; ++r) {
    write_block<FMT_PRHLT_HTK, float>(file, cols, m + r * cols);
  }
}

template <>
void write_matrix<FMT_PRHLT_HTK, double>(
    FILE* file, int rows, int cols, const double* m) {
  for (int r = 0; r < rows; ++r) {
    write_block<FMT_PRHLT_HTK, double>(file, cols, m + r * cols);
  }
}
