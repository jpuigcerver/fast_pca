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
// ---- open_file: Open a file with the specified mode
// ------------------------------------------------------------------------
FILE* open_file(const char* fname, const char* mode);

// ------------------------------------------------------------------------
// ---- open_files: Open a list of files with the specified mode. If the
// ---- list is empty, appends the selected standard file with the given
// ---- name. Useful to read/write from/to multiple files.
// ------------------------------------------------------------------------
void open_files(
    const char* mode, const char* stdname, FILE* stdfile, vector<string>* names,
    vector<FILE*>* files);

// -------------------------------------------------------------------
// ---- read_block: read one data row from an ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
int read_block(FILE* file, int n, real_t* x);

// specialization of read_row function
template <> int read_block<true, float>(FILE*, int, float*);
template <> int read_block<false, float>(FILE*, int, float*);
template <> int read_block<true, double>(FILE*, int, double*);
template <> int read_block<false, double>(FILE*, int, double*);

// -------------------------------------------------------------------
// ---- write_row: write one data row to an ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
void write_block(FILE* file, int n, const real_t* m) {
  if (ascii) {
    for (int c = 0; c < n - 1; ++c) {
      fprintf(file, "%.10g ", m[c]);
    }
    fprintf(file, "%.10g\n", m[n - 1]);
  } else {
    fwrite(m, sizeof(real_t), n, file);
  }
}

#endif  // FAST_PCA_FILE_COMMON_H_
