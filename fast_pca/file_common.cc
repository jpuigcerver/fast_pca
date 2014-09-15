#include "fast_pca/file_common.h"

FILE* file_open(const char* fname, const char* mode) {
  FILE* file = fopen(fname, mode);
  if (!file) {
    fprintf(
        stderr, "ERROR: Failed to open file \"%s\" with mode \"%s\"!\n",
        fname, mode);
    exit(1);
  }
  return file;
}

template <>
int read_row<true, float>(FILE* file, int cols, float* x) {
  int c = 0;
  for (; c < cols && fscanf(file, "%f", x + c) == 1; ++c) {}
  return c;
}

template <>
int read_row<true, double>(FILE* file, int cols, double* x) {
  int c = 0;
  for (; c < cols && fscanf(file, "%lf", x + c) == 1; ++c) {}
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
