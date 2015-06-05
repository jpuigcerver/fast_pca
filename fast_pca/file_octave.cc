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

#include "fast_pca/file_octave.h"

#include <sstream>

using std::istringstream;

template <typename T>
bool read_keyword(FILE* file, const string& key, T* val) {
  bool status = false;
  for (char c = fgetc(file); !feof(file) && !ferror(file); c = fgetc(file)) {
    if (c == '%' || c == '#') {
      // skip whitespace and comment characters.
      for (; !feof(file) && !ferror(file) &&
               (c == ' ' || c == '\t' || c == '%' || c == '#');
           c = fgetc(file)) {}
      string buf = "";
      for (; !feof(file) && !ferror(file) && isalpha(c); c = fgetc(file)) {
        buf += c;
      }
      if (buf == key) {
        // skip whitespace
        for (; !feof(file) && !ferror(file) &&
                 (c == ' ' || c == '\t' || c == ':');
             c = fgetc(file)) {}
        // read keyword value
        buf = "";
        for (; !feof(file) && !ferror(file) && isalnum(c); c = fgetc(file)) {
          buf += c;
        }
        if (!feof(file) && !ferror(file) && !isalpha(c)) {
          ungetc(c, file);
        }
        istringstream iss(buf);
        iss >> *val;
        if (iss) {
          status = true;
        }
        // skip the rest of line
        for (; !feof(file) && !ferror(file) && c != '\n' && c != '\r';
             c = fgetc(file)) {}
        if (!feof(file) && !ferror(file)) {
          c = fgetc(file);
          if (!feof(file) && !ferror(file) && c != '\n' && c != '\r') {
            ungetc(c, file);
          }
        }
        return status;
      }
    }
  }
  return false;
}

// virtual
bool MatrixFile_Octave::copy_header_from(const MatrixFile& other) {
  if (other.format() != format_) return false;
  rows_ = other.rows();
  cols_ = other.cols();
  name_ = static_cast<const MatrixFile_Octave*>(&other)->name_;
  return true;
}

// virtual
bool MatrixFile_Octave::read_header() {
  CHECK(file_);
  string type = "";
  if (!read_keyword(file_, "name", &name_)) return false;
  if (!read_keyword(file_, "type", &type) || type != "matrix") return false;
  if (!read_keyword(file_, "rows", &rows_) || rows_ < 0) return false;
  if (!read_keyword(file_, "columns", &cols_) || cols_ < 0) return false;
  return true;
}

// virtual
void MatrixFile_Octave::write_header() const {
  CHECK(file_);
  if (name_ != "") {
    fprintf(file_, "# name: %s\n", name_.c_str());
  }
  fprintf(file_, "# type: matrix\n");
  fprintf(file_, "# rows: %d\n", rows_);
  fprintf(file_, "# columns: %d\n", cols_);
}

// virtual
int MatrixFile_Octave::read_block(int n, float* m) const {
  CHECK(file_);
  int i = 0;
  for (; i < n && fscanf(file_, "%f", m + i) == 1; ++i) { }
  return i;
}

// virtual
int MatrixFile_Octave::read_block(int n, double* m) const {
  CHECK(file_);
  int i = 0;
  for (; i < n && fscanf(file_, "%lf", m + i) == 1; ++i) { }
  return i;
}

// virtual
void MatrixFile_Octave::write_block(int n, const float* m) const {
  CHECK(file_);
  if (n <= 0) return;
  for (int i = 0; i < n - 1; ++i) fprintf(file_, "%.16g ", m[i]);
  fprintf(file_, "%.16g\n", m[n - 1]);
}

// virtual
void MatrixFile_Octave::write_block(int n, const double* m) const {
  CHECK(file_);
  if (n <= 0) return;
  for (int i = 0; i < n - 1; ++i) fprintf(file_, "%.16g ", m[i]);
  fprintf(file_, "%.16g\n", m[n - 1]);
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_OCTAVE>() {
  return new MatrixFile_Octave;
}

// static
template <>
MatrixFile* MatrixFile::Create<FMT_OCTAVE>(FILE* file) {
  return new MatrixFile_Octave(file);
}
