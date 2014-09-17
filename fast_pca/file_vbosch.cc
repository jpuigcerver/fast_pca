#include "fast_pca/file.h"

template <> int read_matrix_header<FMT_VBOSCH>(
    FILE* file, string* name, int* rows, int* cols) {
  return (fscanf(file, "%d %d", rows, cols) != 2);
}

template <> void write_matrix_header<FMT_VBOSCH>(
    FILE* file, const string& name, int rows, int cols) {
  fprintf(file, "%d %d\n", rows, cols);
}

template <>
int read_block<FMT_VBOSCH, float>(FILE* file, int n, float* m) {
  int i = 0;
  for (; i < n && fscanf(file, "%f", m + i) == 1; ++i) { }
  return i;
}

template <>
int read_block<FMT_VBOSCH, double>(FILE* file, int n, double* m) {
  int i = 0;
  for (; i < n && fscanf(file, "%lf", m + i) == 1; ++i) { }
  return i;
}

template <>
void write_block<FMT_VBOSCH, float>(FILE* file, int n, const float* m) {
  for (int i = 0; i < n - 1; ++i) fprintf(file, "%.12g ", m[i]);
  fprintf(file, "%.12g\n", m[n - 1]);
}

template <>
void write_block<FMT_VBOSCH, double>(FILE* file, int n, const double* m) {
  for (int i = 0; i < n - 1; ++i) fprintf(file, "%.12g ", m[i]);
  fprintf(file, "%.12g\n", m[n - 1]);
}
