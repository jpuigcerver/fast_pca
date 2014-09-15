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
