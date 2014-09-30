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

#include "fast_pca/file.h"

#include <sstream>
#include <string>

using std::istringstream;
using std::string;

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

template <> int read_matrix_header<FMT_OCTAVE>(
    FILE* file, string* name, int* rows, int* cols) {
  string type;
  return !((name == NULL || read_keyword(file, "name", name)) &&
           (read_keyword(file, "type", &type) && type == "matrix") &&
           read_keyword(file, "rows", rows) &&
           read_keyword(file, "columns", cols));
}

template <> void write_matrix_header<FMT_OCTAVE>(
    FILE* file, const string& name, int rows, int cols) {
  if (name != "") {
    fprintf(file, "# name: %s\n", name.c_str());
  }
  fprintf(file, "# type: matrix\n");
  fprintf(file, "# rows: %d\n", rows);
  fprintf(file, "# columns: %d\n", cols);
}

template <>
int read_block<FMT_OCTAVE, float>(FILE* file, int n, float* m) {
  int i = 0;
  for (; i < n && fscanf(file, "%f", m + i) == 1; ++i) { }
  return i;
}

template <>
int read_block<FMT_OCTAVE, double>(FILE* file, int n, double* m) {
  int i = 0;
  for (; i < n && fscanf(file, "%lf", m + i) == 1; ++i) { }
  return i;
}

template <>
void write_block<FMT_OCTAVE, float>(FILE* file, int n, const float* m) {
  for (int i = 0; i < n - 1; ++i) fprintf(file, "%.16g ", m[i]);
  fprintf(file, "%.16g\n", m[n - 1]);
}

template <>
void write_block<FMT_OCTAVE, double>(FILE* file, int n, const double* m) {
  for (int i = 0; i < n - 1; ++i) fprintf(file, "%.16g ", m[i]);
  fprintf(file, "%.16g\n", m[n - 1]);
}

template <>
void write_matrix<FMT_OCTAVE, float>(
    FILE* file, int rows, int cols, const float* m) {
  for (int r = 0; r < rows; ++r) {
    write_block<FMT_OCTAVE, float>(file, cols, m + r * cols);
  }
}

template <>
void write_matrix<FMT_OCTAVE, double>(
    FILE* file, int rows, int cols, const double* m) {
  for (int r = 0; r < rows; ++r) {
    write_block<FMT_OCTAVE, double>(file, cols, m + r * cols);
  }
}

int octave_read_scalar(FILE* file, string* name, int* v) {
  string type;
  return !((name == NULL || read_keyword(file, "name", name)) &&
           (read_keyword(file, "type", &type) && type == "scalar") &&
           fscanf(file, "%d", v) == 1);
}

void octave_write_scalar(FILE* file, const string& name, int n) {
  fprintf(file, "# name: %s\n", name.c_str());
  fprintf(file, "# type: scalar\n");
  fprintf(file, "%d\n", n);
}
