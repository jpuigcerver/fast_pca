#include "file.h"

template <>
int read_row<true, float>(FILE* file, int cols, float* x) {
  int c = 0;
  for (; c < cols && fscanf(file, "%f", x + c) == 1; ++c);
  return c;
}

template <>
int read_row<true, double>(FILE* file, int cols, double* x) {
  int c = 0;
  for (; c < cols && fscanf(file, "%lf", x + c) == 1; ++c);
  return c;
}

template <>
int read_row<false, float>(FILE* file, int cols, float* x) {
  return fread(x, sizeof(float), cols, file);
}

template <>
int read_row<false, double>(FILE* file, int cols, double* x) {
  return fread(x, sizeof(double), cols, file);
}

template <>
void write_row<true, float>(FILE* file, int cols, const float* m) {
  for (int c = 0; c < cols - 1; ++c) {
    fprintf(file, "%g ", m[c]);
  }
  fprintf(file, "%g\n", m[cols - 1]);
}

template <>
void write_row<true, double>(FILE* file, int cols, const double* m) {
  for (int c = 0; c < cols - 1; ++c) {
    fprintf(file, "%lg ", m[c]);
  }
  fprintf(file, "%lg\n", m[cols - 1]);
}

template <>
void write_row<false, float>(FILE* file, int cols, const float* m) {
  fwrite(m, sizeof(float), cols, file);
}

template <>
void write_row<false, double>(FILE* file, int cols, const double* m) {
  fwrite(m, sizeof(double), cols, file);
}


void read_matlab_header(
    char const* fname, FILE* file, int* rows, int* cols) {
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
