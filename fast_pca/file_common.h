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

#ifndef FAST_PCA_FILE_COMMON_H_
#define FAST_PCA_FILE_COMMON_H_

#include <cstdio>
#include <cstdlib>
#include <vector>

using std::vector;

// ------------------------------------------------------------------------
// ---- file_open: Open a file with the specified mode
// ------------------------------------------------------------------------
FILE* file_open(const char* fname, const char* mode);

// -------------------------------------------------------------------
// ---- read_row: read one data row from an ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
int read_row(FILE* file, int cols, real_t* x);

// specialization of read_row function
template <> int read_row<true, float>(FILE*, int, float*);
template <> int read_row<false, float>(FILE*, int, float*);
template <> int read_row<true, double>(FILE*, int, double*);
template <> int read_row<false, double>(FILE*, int, double*);

// -------------------------------------------------------------------
// ---- write_row: write one data row to an ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
void write_row(FILE* file, int cols, const real_t* m) {
  if (ascii) {
    for (int c = 0; c < cols - 1; ++c) {
      fprintf(file, "%.10g ", m[c]);
    }
    fprintf(file, "%.10g\n", m[cols - 1]);
  } else {
    fwrite(m, sizeof(real_t), cols, file);
  }
}


// ------------------------------------------------------------------------
// ---- read_matrix: read ascii/binary matrix from the given file
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void read_matrix(
    const char* fname, FILE* file, int rows, int cols, real_t* m) {
  int tc = 0;
  for (int r = 0; r < rows &&
           (tc = read_row<ascii, real_t>(file, cols, m)) == cols;
       ++r, m += cols) {}
  if (tc != 0 && tc != cols) {
    fprintf(stderr, "ERROR: Corrupted matrix in \"%s\"!\n", fname);
    exit(1);
  }
}


// -------------------------------------------------------------------
// ---- write_matrix: write a matrix to a ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
void write_matrix(FILE* file, int rows, int cols, const real_t* m) {
  for (int r = 0; r < rows; ++r) {
    write_row<ascii, real_t>(file, cols, m + r * cols);
  }
}

#endif  // FAST_PCA_FILE_COMMON_H_