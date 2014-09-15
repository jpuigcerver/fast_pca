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

#ifndef FAST_PCA_FILE_H_
#define FAST_PCA_FILE_H_

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "fast_pca/file_common.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/file_plain.h"
#include "fast_pca/file_octave.h"

using std::vector;

// ------------------------------------------------------------------------
// ---- read_text_header: read MAT header from the given file
// ------------------------------------------------------------------------
void read_text_header(
    char const* fname, FILE* file, int* rows, int* cols);


// ---------------------------------------------------------------------------
// ---- save_text: save a MAT matrix to the given file name in ascii/binary
// ---------------------------------------------------------------------------
template <typename real_t>
void save_text(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = file_open(fname, "w");
  fprintf(file, "%d %d\n", rows, cols);
  write_matrix<true, real_t>(file, rows, cols, m);
  fclose(file);
}


// ------------------------------------------------------------------------
// ---- load_text: load a matrix in MAT format from the given file
// ------------------------------------------------------------------------
template <typename real_t>
void load_text(const char* fname, int* rows, int* cols, vector<real_t>* m) {
  // open MAT file
  FILE* file = file_open(fname, "r");
  read_text_header(fname, file, rows, cols);
  m->resize((*rows) * (*cols));
  read_matrix<true, real_t>(fname, file, rows, *cols, m->data());
  fclose(file);
}




#endif  // FAST_PCA_FILE_H_
