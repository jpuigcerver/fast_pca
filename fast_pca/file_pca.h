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

#include <algorithm>
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
  octave_write_scalar<int>(out_f, "N", n);
  // write m
  write_matrix_header<FMT_OCTAVE>(out_f, "M", 1, d);
  write_block<FMT_OCTAVE, real_t>(out_f, d, m.data());
  // write c
  write_matrix_header<FMT_OCTAVE>(out_f, "C", d, d);
  write_matrix<FMT_OCTAVE, real_t>(out_f, d, d, c.data());
  fclose(out_f);
}

template <typename real_t>
void load_n_mean_cov(
    const string& fname, int* n, int* d, vector<real_t>* m,
    vector<real_t>* c) {
  FILE* file = stdin;
  if (fname != "") { file = open_file(fname.c_str(), "rb"); }
  string ts;
  int tr = -1, tc = -1;
  CHECK_FMT(
      !octave_read_scalar<int>(file, &ts, &tr) && ts == "N",
      "Failed to read scalar N in \"%s\"!", fname.c_str());
  CHECK_FMT(
      *n < 1 || tr == *n,
      "N has a different value (%d) than expected (%d) in file \"%s\"!",
      tr, *n, fname.c_str());
  *n = tr;
  CHECK_FMT(
      !octave_read_matrix(file, &ts, &tr, &tc, m) && ts == "M" && tr == 1,
      "Failed to read matrix M in \"%s\"!", fname.c_str());
  CHECK_FMT(
      *d < 1 || tc == *d,
      "Size of vector M (%d) is different than the expected (%d) in file "
      "\"%s\"!", tc, *d, fname.c_str());
  *d = tc;
  CHECK_FMT(
      !octave_read_matrix(file, &ts, &tr, &tc, c) && ts == "C",
      "Failed to read matrix C in \"%s\"!", fname.c_str());
  CHECK_FMT(
      tr == *d && tc == *d,
      "Size of matrix C (%dx%d) is different than the expected (%dx%d) in "
      "file \"%s\"!", tr, tc, *d, *d, fname.c_str());
  fclose(file);
}

// fname        -> (input) file to store the pca data, "" for stdout
// idim         -> (output) number of dimensions in the original data
// pca_odim     -> (output) number of components available for projection
// exclude_dims -> (output) exclude this number of first/last dimensions
// mean         -> (output) means of each data dimension
// stddev       -> (output) standard deviation of each data dimension
// eigval       -> (output) eigenvalues vector
// eigvec       -> (output) eigenvectors matrix
template <typename real_t>
void load_pca(
    const string& fname, int* idim, int* pca_odim,
    int* exclude_dims, double* remaining_energy, vector<real_t>* mean,
    vector<real_t>* stddev, vector<real_t>* eigval, vector<real_t>* eigvec) {
  FILE* file = stdin;
  if (fname != "") { file = open_file(fname.c_str(), "rb"); }
  string ts;
  int tr = -1, tc = -1;
  // read excluded dimensions
  CHECK_FMT(
      !octave_read_scalar<int>(file, &ts, exclude_dims) && ts == "E",
      "Failed to read E in \"%s\"!", fname.c_str());
  // read remaining energy, not included in the eigenvalues
  CHECK_FMT(
      !octave_read_scalar<double>(file, &ts, remaining_energy) && ts == "R",
      "Failed to read R in \"%s\"!", fname.c_str());
  if (*remaining_energy < 0.0) {
    WARN_FMT("Remaining energy is negative in file \"%s\"...", fname.c_str());
    *remaining_energy = 0.0;
  }
  // read means vector
  CHECK_FMT(
      !octave_read_matrix(file, &ts, &tr, &tc, mean) && ts == "M" && tr == 1,
      "Failed to read M in \"%s\"!", fname.c_str());
  CHECK_FMT(
      *idim < 1 || tc == *idim,
      "Size of vector M (%d) is not the expected (%d) in file \"%s\"!",
      tc, *idim, fname.c_str());
  *idim = tc;
  // read standard deviations vector
  CHECK_FMT(
      !octave_read_matrix(file, &ts, &tr, &tc, stddev) && ts == "S" && tr == 1,
      "Failed to read S in \"%s\"!", fname.c_str());
  CHECK_FMT(
      tc == *idim,
      "Size of vector S (%d) is not the expected (%d) in file \"%s\"!",
      tc, *idim, fname.c_str());
  // read eigenvalues
  CHECK_FMT(
      !octave_read_matrix(file, &ts, &tr, &tc, eigval) && ts == "D" && tr == 1
      && tc > 0, "Failed to read D in \"%s\"!", fname.c_str());
  CHECK_FMT(
      *pca_odim < 0 || *pca_odim <= tc,
      "Size of vector D (%d) is smaller than expected (%d) in file \"%s\"!",
      tc, *pca_odim, fname.c_str());
  if (*pca_odim < 0) *pca_odim = tc;
  // read eigenvectors
  const int pca_idim = *idim - abs(*exclude_dims);
  CHECK_FMT(
      !octave_read_matrix(file, &ts, &tr, &tc, eigvec) && ts == "V",
      "Failed to read V in \"%s\"!", fname.c_str());
  CHECK_FMT(
      tr == pca_idim && *pca_odim <= tc,
      "Size of matrix V (%dx%d) is not the expected (%dx%d) in file \"%s\"!",
      tr, tc, pca_idim, *pca_odim, fname.c_str());
  fclose(file);
}

// fname        -> (input) file to store the pca data, "" for stdout
// exclude_dims -> (input) exclude this number of first/last dimensions
// miss_energy  -> (input) energy not captured by the selected eigenvectors
// mean         -> (input) means of each data dimension
// stddev       -> (input) standard deviation of each data dimension
// eigval       -> (input) eigenvalues vector
// eigvec       -> (input) eigenvectors matrix
template <typename real_t>
void save_pca(
    const string& fname, const int exclude_dims, const double miss_energy,
    const vector<real_t>& mean, const vector<real_t>& stddev,
    const vector<real_t>& eigval, const vector<real_t>& eigvec) {
  // safety checks
  CHECK(mean.size() == stddev.size());
  CHECK(eigval.size() <= mean.size());
  const int pca_idim = mean.size() - abs(exclude_dims);
  const int pca_odim = eigval.size();
  FILE* file = stdout;
  if (fname != "") file = open_file(fname.c_str(), "wb");
  // write excluded dimensions
  octave_write_scalar<int>(file, "E", exclude_dims);
  // remaining energy
  octave_write_scalar<double>(file, "R", miss_energy);
  // write mean
  write_matrix_header<FMT_OCTAVE>(file, "M", 1, mean.size());
  write_block<FMT_OCTAVE, real_t>(file, mean.size(), mean.data());
  // write standard deviation
  write_matrix_header<FMT_OCTAVE>(file, "S", 1, stddev.size());
  write_block<FMT_OCTAVE, real_t>(file, stddev.size(), stddev.data());
  // write eigenvalues
  write_matrix_header<FMT_OCTAVE>(file, "D", 1, eigval.size());
  write_block<FMT_OCTAVE, real_t>(file, eigval.size(), eigval.data());
  // write eigenvectors
  write_matrix_header<FMT_OCTAVE>(file, "V", pca_odim, pca_idim);
  write_matrix<FMT_OCTAVE, real_t>(
      file, pca_idim, pca_odim, eigvec.data());
  fclose(file);
}

#endif  // FAST_PCA_FILE_PCA_H_
