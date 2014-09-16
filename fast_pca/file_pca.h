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

#include "fast_pca/file_common.h"
#include "fast_pca/file_octave.h"

using std::string;
using std::vector;

template <typename real_t>
void load_pca(
    const char* fname, int* d, vector<real_t>* mean, vector<real_t>* stddev,
    vector<real_t>* eigval, vector<real_t>* eigvec) {
  int one = 1;
  FILE* file = file_open(fname, "r");
  // means
  string name = "M";
  octave::read_matrix_header_ascii(fname, file, &name, &one, d);
  mean->resize(*d);
  read_matrix<true, real_t>(fname, file, one, *d, mean->data());
  // standard deviations
  name = "S";
  octave::read_matrix_header_ascii(fname, file, &name, &one, d);
  stddev->resize(*d);
  read_matrix<true, real_t>(fname, file, one, *d, stddev->data());
  // eigenvalues
  name = "D";
  octave::read_matrix_header_ascii(fname, file, &name, &one, d);
  eigval->resize(*d);
  read_matrix<true, real_t>(fname, file, one, *d, eigval->data());
  // eigenvectors
  name = "V";
  octave::read_matrix_header_ascii(fname, file, &name, d, d);
  eigvec->resize(*d);
  read_matrix<true, real_t>(fname, file, *d, *d, eigvec->data());
  fclose(file);
}

template <typename real_t>
void save_pca(
    const char* fname, int d, const vector<real_t>& mean,
    const vector<real_t>& stddev, const vector<real_t>& eigval,
    const vector<real_t>& eigvec) {
  FILE* file = stdout;
  if (fname != NULL) { file = file_open(fname, "w"); }
  // means
  octave::write_matrix_header_ascii(fname, file, "M", 1, d);
  write_row<true, real_t>(file, d, mean.data());
  // standard deviations
  octave::write_matrix_header_ascii(fname, file, "S", 1, d);
  write_row<true, real_t>(file, d, stddev.data());
  // eigenvalues
  octave::write_matrix_header_ascii(fname, file, "D", 1, d);
  write_row<true, real_t>(file, d, eigval.data());
  // eigenvectors
  octave::write_matrix_header_ascii(fname, file, "V", d, d);
  write_matrix<true, real_t>(file, d, d, eigvec.data());
  if (fname != NULL) { fclose(file); }
}

#endif  // FAST_PCA_FILE_PCA_H_
