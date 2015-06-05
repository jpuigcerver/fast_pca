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

#ifndef FAST_PCA_FILE_H_
#define FAST_PCA_FILE_H_

#include <cstdio>
#include <string>
#include <vector>

#include "fast_pca/logging.h"

using std::string;
using std::vector;

typedef enum {
  FMT_UNKNOWN = -1,
  FMT_ASCII   = 0,
  FMT_BINARY  = 1,
  FMT_OCTAVE  = 2,
  FMT_VBOSCH  = 3,
  FMT_HTK     = 4,
  FMT_MAT4    = 5
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


// ------------------------------------------------------------------------
// ---- Abstract templated methods for reading / writing matrices in
// ---- different file formats.
// ---- MatrixHeader class can be inherited in order to add additional
// ---- fields that other formats may need.
// ------------------------------------------------------------------------

class MatrixFile {
 protected:
  FORMAT_CODE format_;
  FILE* file_;
  int rows_;
  int cols_;

 public:
  explicit MatrixFile(FORMAT_CODE format) :
      format_(format), file_(nullptr), rows_(-1), cols_(-1) {}
  explicit MatrixFile(FILE* file) :
      file_(file), rows_(-1), cols_(-1) {}
  virtual ~MatrixFile() {}

  inline FORMAT_CODE format() const { return format_; }
  inline void file(FILE* f) { file_ = f; }
  inline void rows(int n) { rows_ = n; }
  inline void cols(int n) { cols_ = n; }
  inline FILE* file() const { return file_; }
  inline int rows() const { return rows_; }
  inline int cols() const { return cols_; }

  virtual bool read_header() { return true; }
  virtual void write_header() const {}
  virtual bool copy_header_from(const MatrixFile& other) {
    rows_ = other.rows();
    cols_ = other.cols();
    return true;
  }

  virtual int read_block(int n, float* m) const = 0;
  virtual int read_block(int n, double* m) const = 0;
  virtual void write_block(int n, const float* m) const = 0;
  virtual void write_block(int n, const double* m) const = 0;


  template <FORMAT_CODE format>
  static MatrixFile* Create();
  template <FORMAT_CODE format>
  static MatrixFile* Create(FILE* file);
};

#endif  // FAST_PCA_FILE_H_
