#ifndef FILE_H_
#define FILE_H_

#include <cstdio>
#include <cstdlib>
#include <vector>

using std::vector;

// -------------------------------------------------------------------
// ---- read_row: read one data sample from a ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
int read_row(FILE* file, int cols, real_t* x);

// specialization of read_row function
template <> int read_row<true, float>(FILE*, int, float*);
template <> int read_row<false, float>(FILE*, int, float*);
template <> int read_row<true, double>(FILE*, int, double*);
template <> int read_row<false, double>(FILE*, int, double*);

// -------------------------------------------------------------------
// ---- dump_matrix: dump a matrix to a ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
void dump_matrix(FILE* file, int rows, int cols, const real_t* m);

// specialization of dump_matrix function
template <> void dump_matrix<true, float>(FILE*, int, int, const float*);
template <> void dump_matrix<false, float>(FILE*, int, int, const float*);
template <> void dump_matrix<true, double>(FILE*, int, int, const double*);
template <> void dump_matrix<false, double>(FILE*, int, int, const double*);

// ------------------------------------------------------------------------
// ---- read_matlab_header: read MAT header from the given file
// ------------------------------------------------------------------------
void read_matlab_header(const char* fname, FILE* file, int* rows, int* cols) {
  int trows = -1, tcols = -1;
  if (fscanf(file, "%d %d", &trows, &tcols) != 2) {
    fprintf(stderr, "ERROR: Bad MAT header in file \"%s\"!\n", fname);
    exit(1);
  }
  // check matrix size
  if (trows < 1 || tcols < 1) {
    fprintf(stderr, "ERROR: Invalid MAT size in file \"%s\"!\n", fname);
    exit(1);
  }
  if ((*rows > 0 && *rows != trows) || (*cols > 0 && *cols != tcols)) {
    fprintf(stderr, "ERROR: Wrong expected MAT size in file \"%s\"!\n", fname);
    exit(1);
  } else {
    *rows = trows;
    *cols = tcols;
  }
}


// ------------------------------------------------------------------------
// ---- read_matrix: read ascii/binary matrix from the given file
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void read_matrix(
    const char* fname, FILE* file, int* rows, int cols, real_t* m) {
  int tr = 0, tc = 0;
  for (; (tc = read_row<ascii, real_t>(file, cols, m)) == cols;
       ++tr, m += cols);
  if (tc != 0 || tc != cols) {
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
  FILE* file = fopen(fname, ascii ? "w" : "wb");
  if (!file) {
    fprintf(stderr, "ERROR: Failed writing to file \"%s\"!", fname);
    exit(1);
  }
  dump_matrix<ascii, real_t>(file, rows, cols, m);
  fclose(file);
}

// ---------------------------------------------------------------------------
// ---- save_matlab: save a MAT matrix to the given file name in ascii/binary
// ---------------------------------------------------------------------------
template <typename real_t>
void save_matlab(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = fopen(fname, "w");
  if (!file || fprintf(file, "%d %d\n", rows, cols) != 2) {
    fprintf(stderr, "ERROR: Failed writing to file \"%s\"!", fname);
    exit(1);
  }
  dump_matrix<true, real_t>(file, rows, cols, m);
  fclose(file);
}


// ------------------------------------------------------------------------
// ---- load_matrix: load a matrix from the given file name in ascii/binary
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void load_matrix(const char* fname, int rows, int cols, real_t* m) {
  FILE* file = fopen(fname, ascii ? "r" : "rb");
  if (!file) {
    fprintf(stderr, "ERROR: Failed reading file \"%s\"!\n", fname);
    exit(1);
  }
  read_matrix(fname, file, &rows, cols, m);
  fclose(file);
}


// ---------------------------------------------------------------------------
// ---- save_integers: save integer to the given file name in ascii/binary
// ---------------------------------------------------------------------------
template <bool ascii>
void save_integers(const char* fname, int n, const int* v) {
  FILE* file = fopen(fname, ascii ? "w" : "wb");
  if (!file) {
    fprintf(stderr, "ERROR: Failed writing file \"%s\"!\n", fname);
    exit(1);
  }
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
  FILE* file = fopen(fname, ascii ? "r" : "rb");
  if (!file) {
    fprintf(stderr, "ERROR: Failed reading file \"%s\"!\n", fname);
    exit(1);
  }
  if (ascii) {
    for (int i = 0; i < n; ++i) {
      if (fscanf(file, "%d", v + i) != 1) {
        fprintf(stderr, "ERROR: Failed reading integers from \"%s\"!\n", fname);
        exit(1);
      }
    }
  } else {
    if (fread(v, sizeof(int), n, file) != n) {
      fprintf(stderr, "ERROR: Failed reading integers from \"%s\"!\n", fname);
      exit(1);
    }
  }
  fclose(file);
}



// ------------------------------------------------------------------------
// ---- load_matlab: load a matrix in MAT format from the given file
// ------------------------------------------------------------------------
template <typename real_t>
void load_matlab(const char* fname, int* rows, int* cols, vector<real_t>* m) {
  // open MAT file
  FILE* file = fopen(fname, "r");
  if (!file) {
    fprintf(stderr, "ERROR: Failed reading file \"%s\"!\n", fname);
    exit(1);
  }
  read_matlab_header(fname, file, rows, cols);
  m->resize((*rows) * (*cols));
  read_matrix<true, real_t>(fname, file, &rows, cols, m->data());
  fclose(file);
}


#endif  // FILE_H_
