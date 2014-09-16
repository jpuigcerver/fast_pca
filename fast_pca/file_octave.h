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

#ifndef FAST_PCA_FILE_OCTAVE_H_
#define FAST_PCA_FILE_OCTAVE_H_

#include <cstdio>
#include <string>
#include <vector>

#include "fast_pca/file_common.h"

using std::string;
using std::vector;

namespace octave {

void read_matrix_header_ascii(
    const char* fname, FILE* file, string* name, int* rows, int* cols);
void write_matrix_header_ascii(
    const char* fname, FILE* file, const string& name, int rows, int cols);

void read_scalar_ascii(const char* fname, FILE* file, string* name, int* v);
void write_scalar_ascii(
    const char* fname, FILE* file, const string& name, int n);

template <typename real_t>
void save_ascii(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = file_open(fname, "w");
  write_matrix_header_ascii(fname, file, "", rows, cols);
  write_matrix<true, real_t>(file, rows, cols, m);
  fclose(file);
}

template <typename real_t>
void load_ascii(
    const char* fname, int* rows, int* cols, vector<real_t>* m) {
  FILE* file = file_open(fname, "r");
  string name;
  read_matrix_header_ascii(fname, file, &name, rows, cols);
  m->resize((*rows) * (*cols));
  read_matrix<true, real_t>(fname, file, rows, *cols, m->data());
  fclose(file);
}

}  // namespace octave


#endif  // FAST_PCA_FILE_OCTAVE_H_
