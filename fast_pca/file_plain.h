#ifndef FAST_PCA_FILE_PLAIN_H_
#define FAST_PCA_FILE_PLAIN_H_

namespace plain {

template <bool ascii, typename real_t>
void load_plain(const char* fname, int* rows, int cols, real_t* m) {
}

template <bool ascii, typename real_t>
void save_plain(const char* fname, int rows, int cols, real_t* m) {
}

}  // namespace plain

#endif  // FAST_PCA_FILE_PLAIN_H_
