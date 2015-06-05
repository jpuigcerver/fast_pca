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

#ifndef FAST_PCA_FILE_PCA_H_
#define FAST_PCA_FILE_PCA_H_

#include <algorithm>
#include <string>
#include <vector>

#include "fast_pca/file_mat4.h"

template <typename real_t>
void save_n_mean_cov(
    const string& fname, int n, int d, const vector<real_t>& m,
    const vector<real_t>& c) {
  FILE* out_f = stdout;
  if (fname != "") { out_f = open_file(fname.c_str(), "w"); }
  MatrixFile_MAT4::save(out_f, "N", n);
  MatrixFile_MAT4::save(out_f, "M", 1, d, m);
  MatrixFile_MAT4::save(out_f, "C", d, d, c);
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
  int32_t si = 0;
  MatrixFile_MAT4::load(file, &ts, &si);
  CHECK_FMT(
      ts == "N", "Failed to read scalar N in file \"%s\"!", fname.c_str());
  CHECK_FMT(
      *n < 1 || si == *n,
      "N has a different value (%d) than expected (%d) in file \"%s\"!",
      tr, *n, fname.c_str());
  *n = si;
  MatrixFile_MAT4::load(file, &ts, &tr, &tc, m);
  CHECK_FMT(
      ts == "M",
      "Failed to read matrix M in file \"%s\"!", fname.c_str());
  CHECK_FMT(
      tr == 1 && (*d < 1 || tc == *d),
      "Size of vector M (%dx%d) is different than the expected (%dx%d) in file "
      "\"%s\"!", tr, tc, 1, *d, fname.c_str());
  *d = tc;
  MatrixFile_MAT4::load(file, &ts, &tr, &tc, c);
  CHECK_FMT(
      ts == "C",
      "Failed to read matrix C in file \"%s\"!", fname.c_str());
  CHECK_FMT(
      tr == *d && tc == *d,
      "Size of matrix C (%dx%d) is different than the expected (%dx%d) in "
      "file \"%s\"!", tr, tc, *d, *d, fname.c_str());
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
  const int idim = mean.size();
  const int pca_idim = mean.size() - abs(exclude_dims);
  const int pca_odim = eigval.size();
  FILE* file = stdout;
  if (fname != "") file = open_file(fname.c_str(), "wb");
  MatrixFile_MAT4::save(file, "E", static_cast<int32_t>(exclude_dims));
  MatrixFile_MAT4::save(file, "R", miss_energy);
  MatrixFile_MAT4::save(file, "M", idim, 1, mean);
  MatrixFile_MAT4::save(file, "S", idim, 1, stddev);
  MatrixFile_MAT4::save(file, "D", 1, pca_odim, eigval);
  MatrixFile_MAT4::save(file, "V", pca_odim, pca_idim, eigvec);
  fclose(file);
}

// fname        -> (input)  file to store the pca data, "" for stdout
// idim         -> (output) number of dimensions in the original data
// pca_odim     -> (output) number of components available for projection
// exclude_dims -> (output) exclude this number of first/last dimensions
// mean         -> (output) means of each data dimension
// stddev       -> (output) standard deviation of each data dimension
// eigval       -> (output) eigenvalues vector
// eigvec       -> (output) eigenvectors matrix
template <typename real_t>
void load_pca(
    const string& fname, int* exclude_dims, double* remaining_energy,
    vector<real_t>* mean, vector<real_t>* stddev, vector<real_t>* eigval,
    vector<real_t>* eigvec) {
  FILE* file = stdin;
  if (fname != "") { file = open_file(fname.c_str(), "rb"); }
  string ts;
  int tr = -1, tc = -1;
  int32_t si = 0;
  // read excluded dimensions
  MatrixFile_MAT4::load(file, &ts, &si);
  CHECK_FMT(
      ts == "E", "Failed to read E in file \"%s\"!", fname.c_str());
  // read missing energy, not included in the eigenvalues
  MatrixFile_MAT4::load(file, &ts, remaining_energy);
  CHECK_FMT(
      ts == "R", "Failed to read R in file \"%s\"!", fname.c_str());
  if (*remaining_energy < 0.0) {
    WARN_FMT("Remaining energy is negative in file \"%s\"...", fname.c_str());
    *remaining_energy = 0.0;
  }
  // read means vector
  MatrixFile_MAT4::load(file, &ts, &tr, &tc, mean);
  CHECK_FMT(
      ts == "M" && tc == 1 && tr > 0,
      "Failed to read vector M in file \"%s\"!", fname.c_str());
  // read standard deviations vector
  MatrixFile_MAT4::load(file, &ts, &tr, &tc, stddev);
  CHECK_FMT(
      ts == "S" && tc == 1 && tr > 0,
      "Failed to read vector S in file \"%s\"!", fname.c_str());
  CHECK_FMT(
      mean->size() == stddev->size(),
      "Size of vector S is not the same as M in file \"%s\"!", fname.c_str());
  // read eigenvalues
  MatrixFile_MAT4::load(file, &ts, &tr, &tc, eigval);
  CHECK_FMT(
      ts == "D" && tr == 1 && tc >= 0,
      "Failed to read vector D in file \"%s\"!", fname.c_str());
  // read eigenvectors
  MatrixFile_MAT4::load(file, &ts, &tr, &tc, eigvec);
  CHECK_FMT(
      ts == "V" && tr >= 0 && tc >= 0,
      "Failed to read matrix V in file \"%s\"!", fname.c_str());
  fclose(file);
}

#endif  // FAST_PCA_FILE_PCA_H_
