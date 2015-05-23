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

#ifndef FAST_PCA_FAST_PCA_H_
#define FAST_PCA_FAST_PCA_H_

#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/math.h"

using std::string;
using std::vector;

template <FORMAT_CODE fmt, typename real_t>
void compute_mean_comoments_from_inputs(
    int block, vector<string> input, int* n, int* dims,
    vector<real_t>* M, vector<real_t>* C) {
  // open input files
  vector<FILE*> files;
  open_files("rb", "**stdin**", stdin, &input, &files);
  // process files
  const vector<real_t> ones(block, 1);
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int h_rows = -1, h_dims = -1;
    CHECK_FMT(
        !read_matrix_header<fmt>(file, NULL, &h_rows, &h_dims) &&
        (*dims < 0 || h_dims != *dims),
        "Wrong header in matrix file \"%s\"!", fname);
    if (*dims < 0) { *dims = h_dims; }
  }
  CHECK_MSG(*dims > 0, "Number of dimensions couldn't be determined!");
  vector<real_t> x(block * (*dims), 0);  // data block
  vector<real_t> m(*dims, 0);            // mean of the current block
  vector<real_t> d(*dims, 0);            // diff between global and block mean
  M->resize(*dims, 0);                   // resize global mean
  C->resize((*dims) * (*dims), 0);       // resize global co-moments
  *n = 0;                                // total processed rows
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int fr = 0, be = 0, br = 0;
    while ((be = read_block<fmt, real_t>(file, block * *dims, x.data())) > 0) {
      CHECK_FMT(be % (*dims) == 0, "Corrupted matrix in file \"%s\"!\n", fname);
      br = be / (*dims);
      fr += br;
      // compute block mean
      gemv<real_t>(
          'T', br, *dims, 1.0 / br, x.data(), *dims, ones.data(), 1,
          0, m.data(), 1);
      // subtract mean to the current block
      for (int i = 0; i < br; ++i) {
        axpy<real_t>(*dims, -1, m.data(), x.data() + i * (*dims));
      }
      // d = M - m
      memcpy(d.data(), M->data(), sizeof(real_t) * (*dims));
      axpy<real_t>(*dims, -1, m.data(), d.data());
      // update co-moments matrix
      // C += (x - m)' * (x - m)
      gemm<real_t>(
          'T', 'N', *dims, *dims, br, 1, x.data(), *dims, x.data(), *dims, 1,
          C->data(), *dims);
      // C += D * D' * (br * n) / (br + n)
      const int nn = *n + br;
      const real_t cf = br * ((*n) / (1.0 * nn));
      ger<real_t>(*dims, *dims, cf, d.data(), d.data(), C->data());
      // update mean
      for (int i = 0; i < *dims; ++i) {
        (*M)[i] = ((*n) * (*M)[i] + br * m[i]) / nn;
      }
      // update total number of processed rows
      *n = nn;
    }
  }
  // close input files
  close_files(files);
}

template <FORMAT_CODE fmt, typename real_t>
void compute_mean_comoments_from_inputs2(
    int block, vector<string> input, int* n, int* dims,
    vector<real_t>* M, vector<real_t>* C) {
  // open input files
  vector<FILE*> files;
  open_files("rb", "**stdin**", stdin, &input, &files);
  // process files
  const vector<real_t> ones(block, 1);
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int h_rows = -1, h_dims = -1;
    CHECK_FMT(
        !read_matrix_header<fmt>(file, NULL, &h_rows, &h_dims) &&
        (*dims < 0 || h_dims != *dims),
        "Wrong header in matrix file \"%s\"!", fname);
    if (*dims < 0) { *dims = h_dims; }
  }
  CHECK_MSG(*dims > 0, "Number of dimensions couldn't be determined!");
  vector<real_t> x(block * (*dims), 0);  // data block
  vector<real_t> m(*dims, 0);            // mean of the current block
  vector<real_t> d(*dims, 0);            // diff between global and block mean
  M->resize(*dims, 0);                   // resize global mean
  C->resize((*dims) * (*dims), 0);       // resize global co-moments
  *n = 0;                                // total processed rows
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int fr = 0, be = 0, br = 0;
    while ((be = read_block<fmt, real_t>(file, block * *dims, x.data())) > 0) {
      CHECK_FMT(be % (*dims) == 0, "Corrupted matrix in file \"%s\"!\n", fname);
      br = be / (*dims);
      fr += br;
      // compute block mean
      gemv<real_t>(
          'T', br, *dims, 1.0 / br, x.data(), *dims, ones.data(), 1,
          0, m.data(), 1);
      // subtract mean to the current block
      for (int i = 0; i < br; ++i) {
        axpy<real_t>(*dims, -1, m.data(), x.data() + i * (*dims));
      }
      // d = M - m
      memcpy(d.data(), M->data(), sizeof(real_t) * (*dims));
      axpy<real_t>(*dims, -1, m.data(), d.data());
      // update co-moments matrix
      // C += (x - m)' * (x - m)
      gemm<real_t>(
          'T', 'N', *dims, *dims, br, 1, x.data(), *dims, x.data(), *dims, 1.0,
          C->data(), *dims);
      // C += D * D' * (br * n) / (br + n)
      ger<real_t>(
          *dims, *dims, (1.0 * br) * (*n) / ((*n) + br), d.data(), d.data(),
          C->data());
      // update mean
      for (int i = 0; i < *dims; ++i) {
        (*M)[i] = ((*n) * (*M)[i] + br * m[i]) / ((*n) + br);
      }
      // update total number of processed rows
      *n += br;
    }
  }
  // close input files
  close_files(files);
}

#endif  // FAST_PCA_FAST_PCA_H_
