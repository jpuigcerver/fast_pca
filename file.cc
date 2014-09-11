#include "file.h"

template <>
int read_sample<true, float>(FILE* file, int cols, float* x) {
  int c = 0;
  for (; c < cols && fscanf(file, "%f", x + c) == 1; ++c);
  return c;
}

template <>
int read_sample<true, double>(FILE* file, int cols, double* x) {
  int c = 0;
  for (; c < cols && fscanf(file, "%lf", x + c) == 1; ++c);
  return c;
}

template <>
int read_sample<false, float>(FILE* file, int cols, float* x) {
  return fread(x, sizeof(float), cols, file);
}

template <>
int read_sample<false, double>(FILE* file, int cols, double* x) {
  return fread(x, sizeof(double), cols, file);
}

template <>
void dump_matrix<true, float>(
    FILE* file, int rows, int cols, const float* m) {
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      fprintf(file, "%g ", m[r * cols + c]);
    }
    fprintf(file, "\n");
  }
}

template <>
void dump_matrix<true, double>(
    FILE* file, int rows, int cols, const double* m) {
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) {
      fprintf(file, "%lg ", m[r * cols + c]);
    }
    fprintf(file, "\n");
  }
}

template <>
void dump_matrix<false, float>(
    FILE* file, int rows, int cols, const float* m) {
  fwrite(m, sizeof(float), rows * cols, file);
}

template <>
void dump_matrix<false, double>(
    FILE* file, int rows, int cols, const double* m) {
  fwrite(m, sizeof(double), rows * cols, file);
}
