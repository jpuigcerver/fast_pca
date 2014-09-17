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
#include <string>
#include <vector>

using std::string;
using std::vector;

typedef enum {
  FMT_UNKNOWN = -1,
  FMT_ASCII = 0,
  FMT_BINARY,
  FMT_OCTAVE,
  FMT_VBOSCH
} FORMAT_CODE;

FORMAT_CODE format_code_from_name(const string& name);

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

// ------------------------------------------------------------------------
// ---- close_files: Close a list of opened files
// ------------------------------------------------------------------------
void close_files(const vector<FILE*>& files);

template <int code>
int read_matrix_header(FILE* file, string* name, int* rows, int* cols);

template <int code>
void write_matrix_header(FILE* file, const string& name, int rows, int cols);

template <int code, typename real_t>
int read_block(FILE* file, int n, real_t* m);

template <int code, typename real_t>
void write_block(FILE* file, int n, const real_t* m);


#endif  // FAST_PCA_FILE_H_
