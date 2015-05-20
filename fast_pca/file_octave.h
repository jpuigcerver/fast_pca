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

#include <string>
#include <vector>

using std::string;
using std::vector;

int octave_read_int(FILE* file, string* name, int* v);
void octave_write_int(FILE* file, const string& name, int n);

template <typename real_t>
int octave_read_matrix(
    FILE* file, string* name, int* r, int* c, vector<real_t>* m) {
  if (read_matrix_header<FMT_OCTAVE>(file, name, r, c) ||
      *r <= 0 || *c <= 0) return 1;
  const int s = (*r) * (*c);
  m->resize(s);
  return (read_block<FMT_OCTAVE, real_t>(file, s, m->data()) != s);
}

#endif  // FAST_PCA_FILE_OCTAVE_H_
