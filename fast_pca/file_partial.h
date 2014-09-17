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

#ifndef FAST_PCA_FILE_PARTIAL_H_
#define FAST_PCA_FILE_PARTIAL_H_

#include <string>
#include <vector>

#include "fast_pca/file_common.h"
#include "fast_pca/file_octave.h"

using std::string;
using std::vector;

namespace partial {

template <typename real_t>
void save_partial(
    const char* fname, int n, int d, const vector<real_t>& m,
    const vector<real_t>& c) {
  FILE* file = open_file(fname, "w");
  /*// write number of data rows processed
  octave::write_scalar_ascii(fname, file, "N", n);
  // write partial result for the mean
  octave::write_matrix_header_ascii(fname, file, "M", 1, d);
  write_row<true, real_t>(file, d, m.data());
  // write partial result for the covariance matrix
  octave::write_matrix_header_ascii(fname, file, "C", d, d);
  write_matrix<true, real_t>(file, d, d, c.data());*/
  fclose(file);
}

template <typename real_t>
void load_partial(
    const char* fname, int* n, int* d, vector<real_t>* m, vector<real_t>* c) {
  int one = 1;
  FILE* file = open_file(fname, "r");
  /*
  // read number of processed rows
  string name = "N";
  octave::read_scalar_ascii(fname, file, &name, n);
  // read partial result for the mean
  name = "M";
  octave::read_matrix_header_ascii(fname, file, &name, &one, d);
  m->resize(*d);
  read_row<true, real_t>(file, *d, m->data());
  // read partial result for the covariance
  name = "C";
  octave::read_matrix_header_ascii(fname, file, &name, d, d);
  c->resize((*d) * (*d));
  read_matrix<true, real_t>(fname, file, *d, *d, c->data());
  */
  fclose(file);
}

}  // namespace partial

#endif  // FAST_PCA_FILE_PARTIAL_H_
