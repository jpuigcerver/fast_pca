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

#include "fast_pca/file_common.h"

FILE* open_file(const char* fname, const char* mode) {
  FILE* file = fopen(fname, mode);
  if (!file) {
    fprintf(
        stderr, "ERROR: Failed to open file \"%s\" with mode \"%s\"!\n",
        fname, mode);
    exit(1);
  }
  return file;
}

void open_files(
    const char* mode, const char* stdname, FILE* stdfile, vector<string>* names,
    vector<FILE*>* files) {
  if (names->empty()) {
    names->push_back(stdname);
    files->push_back(stdfile);
  } else {
    for (size_t f = 0; f < names->size(); ++f) {
      files->push_back(open_file((*names)[f].c_str(), mode));
    }
  }
}

void close_files(const vector<FILE*>& files) {
  for (size_t f = 0; f < files.size(); ++f) {
    FILE* file = files[f];
    if (file != stdin && file != stdout && file != stderr) {
      fclose(file);
    }
  }
}

template <>
int read_block<true, float>(FILE* file, int n, float* x) {
  int c = 0;
  for (; c < n && fscanf(file, "%f", x + c) == 1; ++c) {}
  return c;
}

template <>
int read_block<true, double>(FILE* file, int n, double* x) {
  int c = 0;
  for (; c < n && fscanf(file, "%lf", x + c) == 1; ++c) {}
  return c;
}

template <>
int read_block<false, float>(FILE* file, int n, float* x) {
  return fread(x, sizeof(float), n, file);
}

template <>
int read_block<false, double>(FILE* file, int n, double* x) {
  return fread(x, sizeof(double), n, file);
}
