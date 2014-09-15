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
#include <cstdlib>
#include <vector>

using std::vector;

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
      fprintf(file, "%.8g ", m[c]);
    }
    fprintf(file, "%.8g\n", m[cols - 1]);
  } else {
    fwrite(m, sizeof(real_t), cols, file);
  }
}


// -------------------------------------------------------------------
// ---- dump_matrix: dump a matrix to a ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
void dump_matrix(FILE* file, int rows, int cols, const real_t* m) {
  for (int r = 0; r < rows; ++r) {
    write_row<ascii, real_t>(file, cols, m + r * cols);
  }
}


// ------------------------------------------------------------------------
// ---- file_open: Open a file with the specified mode
// ------------------------------------------------------------------------
FILE* file_open(const char* fname, const char* mode);


// ------------------------------------------------------------------------
// ---- read_text_header: read MAT header from the given file
// ------------------------------------------------------------------------
void read_text_header(
    char const* fname, FILE* file, int* rows, int* cols);


// ------------------------------------------------------------------------
// ---- read_matrix: read ascii/binary matrix from the given file
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void read_matrix(
    const char* fname, FILE* file, int* rows, int cols, real_t* m) {
  int tr = 0, tc = 0;
  for (; (tc = read_row<ascii, real_t>(file, cols, m)) == cols;
       ++tr, m += cols) {}
  if (tc != 0 && tc != cols) {
    fprintf(stderr, "ERROR: Corrupted matrix in \"%s\"!\n", fname);
    exit(1);
  }
  if (*rows > 0 && tr != *rows) {
    fprintf(stderr, "ERROR: Matrix in \"%s\" is too big!\n", fname);
    exit(1);
  } else {
    *rows = tr;
  }
}

// ------------------------------------------------------------------------
// ---- save_matrix: save a matrix to the given file name in ascii/binary
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void save_matrix(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = file_open(fname, ascii ? "w" : "wb");
  dump_matrix<ascii, real_t>(file, rows, cols, m);
  fclose(file);
}


// ---------------------------------------------------------------------------
// ---- save_text: save a MAT matrix to the given file name in ascii/binary
// ---------------------------------------------------------------------------
template <typename real_t>
void save_text(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = file_open(fname, "w");
  fprintf(file, "%d %d\n", rows, cols);
  dump_matrix<true, real_t>(file, rows, cols, m);
  fclose(file);
}


// ------------------------------------------------------------------------
// ---- load_matrix: load a matrix from the given file name in ascii/binary
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void load_matrix(const char* fname, int rows, int cols, real_t* m) {
  FILE* file = file_open(fname, ascii ? "w" : "wb");
  read_matrix(fname, file, &rows, cols, m);
  fclose(file);
}


// ---------------------------------------------------------------------------
// ---- save_integers: save integer to the given file name in ascii/binary
// ---------------------------------------------------------------------------
template <bool ascii>
void save_integers(const char* fname, int n, const int* v) {
  FILE* file = file_open(fname, ascii ? "w" : "wb");
  if (ascii) {
    for (int i = 0; i < n; ++i) {
      fprintf(file, "%d ", v[i]);
    }
  } else {
    fwrite(v, sizeof(int), n, file);
  }
  fclose(file);
}


// ---------------------------------------------------------------------------
// ---- load_integers: load integers from the given file name in ascii/binary
// ---------------------------------------------------------------------------
template <bool ascii>
void load_integers(const char* fname, int n, int* v) {
  FILE* file = file_open(fname, ascii ? "r" : "rb");
  if (ascii) {
    for (int i = 0; i < n; ++i) {
      if (fscanf(file, "%d", v + i) != 1) {
        fprintf(stderr, "ERROR: Failed reading integers from \"%s\"!\n", fname);
        exit(1);
      }
    }
  } else {
    if (fread(v, sizeof(int), n, file) != static_cast<size_t>(n)) {
      fprintf(stderr, "ERROR: Failed reading integers from \"%s\"!\n", fname);
      exit(1);
    }
  }
  fclose(file);
}


// ------------------------------------------------------------------------
// ---- load_text: load a matrix in MAT format from the given file
// ------------------------------------------------------------------------
template <typename real_t>
void load_text(const char* fname, int* rows, int* cols, vector<real_t>* m) {
  // open MAT file
  FILE* file = file_open(fname, "r");
  read_text_header(fname, file, rows, cols);
  m->resize((*rows) * (*cols));
  read_matrix<true, real_t>(fname, file, rows, *cols, m->data());
  fclose(file);
}




#endif  // FAST_PCA_FILE_H_
