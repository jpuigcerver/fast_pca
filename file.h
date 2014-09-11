#ifndef FILE_H_
#define FILE_H_

#include <cstdio>
#include <cstdlib>

// -------------------------------------------------------------------
// ---- read_sample: read one data sample from a ascii/binary file
// -------------------------------------------------------------------
template <bool ascii, typename real_t>
int read_sample(FILE* file, int cols, real_t* x);

// specialization of read_sample function
template <> int read_sample<true, float>(FILE*, int, float*);
template <> int read_sample<false, float>(FILE*, int, float*);
template <> int read_sample<true, double>(FILE*, int, double*);
template <> int read_sample<false, double>(FILE*, int, double*);

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
// ---- save_matrix: save a matrix to the given file name in ascii/binary
// ------------------------------------------------------------------------
template <bool ascii, typename real_t>
void save_matrix(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = fopen(fname, ascii ? "w" : "wb");
  if (!file) {
    fprintf(stderr, "ERROR: Failed writing file \"%s\"!", fname);
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
  if (!file) {
    fprintf(stderr, "ERROR: Failed writing file \"%s\"!", fname);
    exit(1);
  }
  fprintf(file, "%d %d\n", rows, cols);
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
    fprintf(stderr, "ERROR: Failed writing file \"%s\"!\n", fname);
    exit(1);
  }
  for (int r = 0; r < rows; ++r) {
    if (read_sample<ascii, real_t>(file, cols, m + r * cols) != cols) {
      fprintf(stderr, "ERROR: Bad matrix format in \"%s\"!\n", fname);
      exit(1);
    }
  }
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
    fprintf(stderr, "ERROR: Failed writing file \"%s\"!\n", fname);
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


#endif  // FILE_H_
