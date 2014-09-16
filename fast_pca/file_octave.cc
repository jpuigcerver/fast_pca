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

#include "fast_pca/file_octave.h"

#include <sstream>
#include <string>

using std::istringstream;
using std::string;

namespace octave {

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

void read_matrix_header_ascii(
    const char* fname, FILE* file, string* name, int* rows, int* cols) {
  int tr = -1, tc = -1;
  string tn, tm;
  if (!read_keyword(file, "name", &tn) ||
      !read_keyword(file, "type", &tm) ||
      !read_keyword(file, "rows", &tr) ||
      !read_keyword(file, "columns", &tc)) {
    fprintf(stderr, "ERROR: Bad Octave ASCII header in \"%s\"!\n", fname);
    exit(1);
  }
  if (tm != "matrix") {
    fprintf(
        stderr, "ERROR: Ivalid matrix type: \"%s\" (expected: \"matrix\")!\n",
        tm.c_str());
    exit(1);
  }
  if (*name != "" && tn != *name) {
    fprintf(
        stderr, "ERROR: Wrong matrix name: \"%s\" (expected: \"%s\")!\n",
        tn.c_str(), name->c_str());
  }
  if ((*rows > 0 && *rows != tr) || (*cols > 0 && *cols != tc)) {
    fprintf(stderr, "ERROR: Wrong matrix size read from \"%s\"!\n", fname);
    exit(1);
  }
  *rows = tr;
  *cols = tc;
  *name = tn;
}

void read_scalar_ascii(
    const char* fname, FILE* file, string* name, int* v) {
  string tn, tm;
  if (!read_keyword(file, "name", &tn) ||
      !read_keyword(file, "type", &tm) ||
      fscanf(file, "%d", v) != 1) {
    fprintf(stderr, "ERROR: Bad Octave ASCII scalar in \"%s\"!\n", fname);
    exit(1);
  }
  if (tm != "scalar") {
    fprintf(stderr, "ERROR: Invalid scalar type: \"%s\"!\n", tm.c_str());
    exit(1);
  }
  if (*name != "" && tn != *name) {
    fprintf(
        stderr, "ERROR: Wrong matrix name: \"%s\" (expected: \"%s\")!\n",
        tn.c_str(), name->c_str());
  }
}

void write_matrix_header_ascii(
    const char* fname, FILE* file, const string& name, int rows, int cols) {
  if (name != "") {
    fprintf(file, "# name: %s\n", name.c_str());
  }
  fprintf(file, "# type: matrix\n");
  fprintf(file, "# rows: %d\n", rows);
  fprintf(file, "# columns: %d\n", cols);
}

void write_scalar_ascii(
    const char* fname, FILE* file, const string& name, int n) {
  if (name != "") {
    fprintf(file, "# name: %s\n", name.c_str());
  }
  fprintf(file, "# type: scalar\n");
  fprintf(file, "%d\n", n);
}

}  // namespace octave
