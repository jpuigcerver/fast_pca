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

#ifndef FAST_PCA_PARTIAL_H_
#define FAST_PCA_PARTIAL_H_

#include <numeric>
#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/math.h"

using std::accumulate;
using std::string;
using std::vector;

// Compute sufficient variables to compute the covariance matrix
template <bool ascii, typename real_t>
int compute_partial(
    FILE* file, int block_size, int dims, real_t* C, real_t* S) {
  const vector<real_t> ones(block_size, 1);
  vector<real_t> x(block_size * dims);
  int d = 0, b = 0, n = 0;
  do {
    const int b = read_block<ascii, real_t>(file, block_size * dims, x.data());
    if (b % dims != 0) {
      fprintf(stderr, "ERROR: Corrupted data matrix!\n");
      exit(1);
    }
    const int read_rows = b / dims;
    // update covariance partial
    gemm<real_t>(
        'T', 'N', dims, read_rows, dims, 1, x.data(), read_rows, x.data(), dims,
        1, C, dims);
    // update mean partial
    gemv<real_t>('T', dims, read_rows, 1, x.data(), read_rows,
  } while(b == block_size);

  for (; (d = read_row<ascii, real_t>(file, dims, x.data())) == dims; ++n) {
    axpy<real_t>(dims, 1.0, x.data(), S);
    ger<real_t>(dims, dims, 1.0, x.data(), x.data(), C);
  }

  return (d == 0 || d == dims) ? n : -1;
}

template <typename real_t>
int reduce_partial(int mappers, int dims, int* n, real_t* C, real_t* S) {
  const int N = accumulate(n, n + mappers, 0);
  for (int m = 1; m < mappers; ++m) {
    axpy<real_t>(dims, 1.0, S + m * dims, S);
    axpy<real_t>(dims * dims, 1.0, C + m * dims * dims, C);
  }
  return N;
}

template <typename real_t>
int partial_cov_mean(
    const string& format, int *dims, const vector<string>& input,
    const string& output, vector<real_t>* m, vector<real_t>* c) {
  const bool ascii = format == "binary" ? false : true;
  // open input files
  vector<FILE*> input_files;
  if (input.size() == 0) { input_files.push_back(stdin); }
  for (size_t f = 0; f < input.size(); ++f) {
    FILE* file = file_open(input[f].c_str(), format == "binary" ? "rb" : "r");
    input_files.push_back(file);
  }
  // number of expected rows in each input file
  vector<int> expect_rows(input_files.size(), -1);
  // number of processed rows in each input file
  vector<int> processed_rows(input_files.size(), 0);
  // read input headers
  if (format == "simple") {
    for (size_t f = 0; f < input_files.size(); ++f) {
      FILE* file = input_files[f];
      const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
      simple::read_header_ascii(fname, file, &expect_rows[f], dims);
    }
  } else if (*dims < 1) {
    fprintf(stderr, "ERROR: You must specify the input dimensions!\n");
    exit(1);
  }
  // this vector accumulates the sum across each dimension, needed to compute
  // the mean and the covariance of the data.
  m->resize(*dims, 0);
  // this matrix accumulates the sum of the cross-product across dimensions,
  // needed to compute the covariance of the data.
  c->resize((*dims) * (*dims), 0);
  // TODO(jpuigcerver): this can run in parallel
  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
    // compute partial results
    processed_rows[f] = ascii ?
        compute_partial<true, real_t>(
            file, *dims, c->data(), m->data()) :
        compute_partial<false, real_t>(
            file, *dims, c->data(), m->data());
    if (processed_rows[f] < 0 ||
        (expect_rows[f] > 0 && processed_rows[f] != expect_rows[f])) {
      fprintf(stderr, "ERROR: Corrupted matrix file \"%s\"!\n", fname);
      exit(1);
    }
  }
  return reduce_partial<real_t>(
      input_files.size(), *dims, processed_rows.data(), c->data(), m->data());
}

#endif  // FAST_PCA_PARTIAL_H_
