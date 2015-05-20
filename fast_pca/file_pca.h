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

#ifndef FAST_PCA_FILE_PCA_H_
#define FAST_PCA_FILE_PCA_H_

#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/file_octave.h"
#include "fast_pca/logging.h"

using std::string;
using std::vector;

template <typename real_t>
void save_n_mean_cov(
    const string& fname, int n, int d, const vector<real_t>& m,
    const vector<real_t>& c) {
  FILE* out_f = stdout;
  if (fname != "") { out_f = open_file(fname.c_str(), "w"); }
  // write n
  octave_write_int(out_f, "N", n);
  // write m
  write_matrix_header<FMT_OCTAVE>(out_f, "M", 1, d);
  write_block<FMT_OCTAVE, real_t>(out_f, d, m.data());
  // write c
  write_matrix_header<FMT_OCTAVE>(out_f, "C", d, d);
  write_matrix<FMT_OCTAVE, real_t>(out_f, d, d, c.data());
  if (fname != "") { fclose(out_f); }
}

template <typename real_t>
void load_n_mean_cov(
    const string& fname, int* n, int* d, vector<real_t>* m,
    vector<real_t>* c) {
  FILE* file = stdin;
  if (fname != "") { file = open_file(fname.c_str(), "rb"); }
  string ts;
  int tr = -1, tc = -1;
  if (octave_read_int(file, &ts, &tr) || ts != "N" ||
      (*n >= 0 && tr != *n)) {
    ERROR("Failed to read scalar N in \"%s\"!", fname.c_str());
  }
  *n = tr;
  if (octave_read_matrix(file, &ts, &tr, &tc, m) || ts != "M" ||
      tr != 1 || (*d >= 0 && tc != *d)) {
    ERROR("Failed to read matrix M in \"%s\"!", fname.c_str());
  }
  *d = tc;
  if (octave_read_matrix(file, &ts, &tr, &tc, c) || ts != "C" ||
      tr != *d || tc != *d) {
    ERROR("Failed to read matrix C in \"%s\"!", fname.c_str());
  }
  if (fname != "") { fclose(file); }
}

template <typename real_t>
void load_pca(
    const string& fname, int* d, int* ed, vector<real_t>* mean,
    vector<real_t>* stddev, vector<real_t>* eigval, vector<real_t>* eigvec) {
  FILE* file = stdin;
  if (fname != "") { file = open_file(fname.c_str(), "rb"); }
  string ts;
  int tr = -1, tc = -1;
  // read excluded dimensions, if it cannot be read, assume 0.
  if (octave_read_int(file, &ts, ed)) *ed = 0;
  // read means vector
  CHECK_MSG(
      !octave_read_matrix(file, &ts, &tr, &tc, mean) && ts == "M" && tr == 1 &&
      (*d < 1 || tc == *d), "Failed to read M in \"%s\"!", fname.c_str());
  *d = tc;
  const int effective_dims = *d - abs(*ed);
  CHECK_MSG(
      effective_dims > 0, "Too many excluded dimensions "
      "(dims = %d, excluded = %d)!", *d, *ed);
  // read standard deviations vector
  CHECK_MSG(
      !octave_read_matrix(file, &ts, &tr, &tc, stddev) && ts == "S" &&
      tr == 1 && tc == *d, "Failed to read S in \"%s\"!", fname.c_str());
  // read eigenvalues
  CHECK_MSG(
      !octave_read_matrix(file, &ts, &tr, &tc, eigval) && ts == "D" &&
      tr == 1 && tc == effective_dims, "Failed to read D in \"%s\"!",
      fname.c_str());
  // read eigenvectors
  CHECK_MSG(
      !octave_read_matrix(file, &ts, &tr, &tc, eigvec) && ts == "V" &&
      tr == effective_dims && tc == effective_dims,
      "Failed to read V in \"%s\"!", fname.c_str());
  if (fname != "") { fclose(file); }
}

template <typename real_t>
void save_pca(
    const string& fname, int d, int exclude_dims, const vector<real_t>& mean,
    const vector<real_t>& stddev, const vector<real_t>& eigval,
    const vector<real_t>& eigvec) {
  FILE* file = stdout;
  if (fname != "") file = open_file(fname.c_str(), "wb");
  // write excluded dimensions
  octave_write_int(file, "E", exclude_dims);
  // dimensions of the projection matrix and the eigenvalues
  const int ed = d - abs(exclude_dims);
  // write mean
  write_matrix_header<FMT_OCTAVE>(file, "M", 1, d);
  write_block<FMT_OCTAVE, real_t>(file, d, mean.data());
  // write standard deviation
  write_matrix_header<FMT_OCTAVE>(file, "S", 1, d);
  write_block<FMT_OCTAVE, real_t>(file, d, stddev.data());
  // write eigenvalues
  write_matrix_header<FMT_OCTAVE>(file, "D", 1, ed);
  write_block<FMT_OCTAVE, real_t>(file, ed, eigval.data());
  // write eigenvectors
  write_matrix_header<FMT_OCTAVE>(file, "V", ed, ed);
  write_matrix<FMT_OCTAVE, real_t>(file, ed, ed, eigvec.data());
  if (fname != "") fclose(file);
}

#endif  // FAST_PCA_FILE_PCA_H_
