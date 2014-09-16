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

#ifndef FAST_PCA_FILE_PLAIN_H_
#define FAST_PCA_FILE_PLAIN_H_

#include "fast_pca/file_common.h"

namespace plain {

template <bool ascii, typename real_t>
void load_plain(const char* fname, int rows, int cols, real_t* m) {
  FILE* file = file_open(fname, ascii ? "r" : "rb");
  read_matrix<ascii, real_t>(fname, file, rows, cols, m);
  fclose(file);
}

template <bool ascii, typename real_t>
void save_plain(const char* fname, int rows, int cols, real_t* m) {
  FILE* file = file_open(fname, ascii ? "w" : "wb");
  write_matrix<ascii, real_t>(fname, file, rows, cols, m);
  fclose(file);
}

}  // namespace plain

#endif  // FAST_PCA_FILE_PLAIN_H_
