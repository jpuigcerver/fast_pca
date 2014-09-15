#include "fast_pca/file_octave.h"

#include <string>

using std::string;

namespace octave {

bool read_keyword_int(FILE *file, const char *key, int* val) {
  bool status = false;
  for(char c = fgetc(file); !feof(file) && !ferror(file); c = fgetc(file)) {
    if (c == '%' || c == '#') {
      // skip whitespace and comment characters.
      for (; !feof(file) && !ferror(file) &&
               (c == ' ' || c == '\t' || c == '%' || c == '#');
           c = fgetc(file)) {}
      string buf;
      for (; !feof(file) && !ferror(file) && isalpha(c); c = fgetc(file)) {
        buf += c;
      }
      if (buf == key) {
        // skip whitespace
        for (; !feof(file) && !ferror(file) && (c == ' ' || c == '\t');
             c = fgetc(file)) {}
        // read keyword value
        if (!feof(file) && !ferror(file) && c != '\n' && c != '\r' &&
            ungetc(c, file) == c && fscanf(file, "%d", val) == 1) {
          status = true;
        }
        // skip the rest of line
        for (; !feof(file) && !ferror(file) && c != '\n' && c != '\r';
             c = fgetc(file)) {}
        if (!feof(file) && !ferror(file)) {
          c = fgetc(file);
          if (!feof(file) && !ferror(file) && c != '\n' && c != '\r') {
            ungetc(c, file);
          }
        }
        return status;
      }
    }
  }
  return false;
}

void read_header_ascii(
    const char* fname, FILE* file, int* rows, int* cols) {
  int tr = -1, tc = -1;
  if (!read_keyword_int(file, "rows:", &tr) ||
      !read_keyword_int(file, "columns:", &tc)) {
    fprintf(stderr, "ERROR: Bad Octave ASCII file in \"%s\"!\n", fname);
    exit(1);
  }
  if ((*rows > 0 && *rows != tr) || (*cols > 0 && *cols != tr)) {
    fprintf(stderr, "ERROR: Wrong matrix size read from \"%s\"!\n", fname);
    exit(1);
  }
  *rows = tr;
  *cols = tc;
}

void write_header_ascii(const char* fname, FILE* file, int rows, int cols) {
  fprintf(file, "# type: matrix\n");
  fprintf(file, "# rows: %d\n", rows);
  fprintf(file, "# columns: %d\n", cols);
}

}  // namespace octave
