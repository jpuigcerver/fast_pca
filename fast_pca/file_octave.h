#ifndef FAST_PCA_FILE_OCTAVE_H_
#define FAST_PCA_FILE_OCTAVE_H_

namespace octave {

void read_header_ascii(const char* fname, FILE* file, int* rows, int* cols);
void write_header_ascii(const char* fname, FILE* file, int rows, int cols);

template <typename real_t>
void save_octave_ascii(const char* fname, int rows, int cols, const real_t* m) {
  FILE* file = file_open(fname, "w");
  write_header_ascii(fname, file, rows, cols);
  dump_matrix<true, real_t>(file, rows, cols, m);
  fclose(file);
}

template <typename real_t>
void load_octave_ascii(
    const char* fname, int* rows, int* cols, vector<real_t>* m) {
  FILE* file = file_open(fname, "r");
  read_header_ascii(fname, file, rows, cols);
  m->resize((*rows) * (*cols));
  read_matrix<true, real_t>(fname, file, rows, *cols, m->data());
  fclose(file);
}

}  // namespace octave


#endif  // FAST_PCA_FILE_OCTAVE_H_
