#ifndef FAST_PCA_FILE_PCA_H_
#define FAST_PCA_FILE_PCA_H_

#include "fast_pca/file_common.h"
#include "fast_pca/file_octave.h"

template <typename real_t>
void load_pca(
    const char* fname, int* d, vector<real_t>* mean, vector<real_t>* stddev,
    vector<real_t>* eigval, vector<real_t>* eigvec) {
  int one = 1;
  FILE* file = file_open(fname, "r");
  // means
  octave::read_header_ascii(fname, file, &one, d);
  mean->reserve(*d);
  read_matrix<true, real_t>(fname, file, one, *d, mean->data());
  // standard deviations
  octave::read_header_ascii(fname, file, &one, d);
  stddev->reserve(*d);
  read_matrix<true, real_t>(fname, file, one, *d, stddev->data());
  // eigenvalues
  octave::read_header_ascii(fname, file, &one, d);
  eigval->reserve(*d);
  read_matrix<true, real_t>(fname, file, one, *d, eigval->data());
  // eigenvectors
  octave::read_header_ascii(fname, file, d, d);
  eigvec->reserve(*d);
  read_matrix<true, real_t>(fname, file, *d, *d, eigvec->data());
  fclose(file);
}

template <typename real_t>
void save_pca(
    const char* fname, int d, const vector<real_t>& mean,
    const vector<real_t>& stddev, const vector<real_t>& eigval,
    const vector<real_t>& eigvec) {
  FILE* file = stdout;
  if (fname != NULL) { file = file_open(fname, "w"); }
  // means
  fprintf(file, "# name: M\n");
  octave::write_header_ascii(fname, file, 1, d);
  write_row<true, real_t>(file, d, mean.data());
  // standard deviations
  fprintf(file, "# name: S\n");
  octave::write_header_ascii(fname, file, 1, d);
  write_row<true, real_t>(file, d, stddev.data());
  // eigenvalues
  fprintf(file, "# name: D\n");
  octave::write_header_ascii(fname, file, 1, d);
  write_row<true, real_t>(file, d, eigval.data());
  // eigenvectors
  fprintf(file, "# name: V\n");
  octave::write_header_ascii(fname, file, d, d);
  write_matrix<true, real_t>(file, d, d, eigvec.data());
  if (fname != NULL) { fclose(file); }
}

#endif  // FAST_PCA_FILE_PCA_H_
