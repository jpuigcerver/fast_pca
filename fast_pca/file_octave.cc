#include "fast_pca/file_octave.h"
#include "fast_pca/file_common.h"

namespace octave {

bool read_keyword_int(FILE *file, const char *key, int* val) {
  bool status = false;
  for(char c = fgetc(file); !feof(file) && !ferror(file); c = fgetc(file)) {
    if (c == '%' || c == '#') {
      // skip whitespace and comment characters.
      while(; !feof(file) && !ferror(file) &&
            (c == ' ' || c == '\t' || c == '%' || c == '#');
            c = fgetc(file)) {}
      string buf;
      while(; !feof(file) && !ferror(file) && isalpha(c); c = fgetc(file)) {
        buf += c;
      }
      if (buf == key) {
        // skip whitespace
        while(; !feof(file) && !ferror(file) && (c == ' ' || c == '\t');
              c = fgetc(file)) {}
        // read keyword value
        if (!feof(file) && !ferror(file) && c != '\n' && c != '\r' &&
            ungetc(c, file) && fscanf("%d", val)) {
          status = true;
        }
        // skip the rest of line
      }
    }
  }
  return status;
}

void read_header_ascii(const char* fname, FILE* file, int* rows, int* cols) {
  if (!read_keyword_int(file, "rows:", rows) ||
      !read_keyword_int(file, "columns:", cols)) {
    fprintf(stderr, "ERROR: Bad Octave ACSII file in \"%s\"!\n", fname);
    exit(1);
  }
}

void write_header_ascii(const char* fname, FILE* file, int rows, int cols) {
  fprintf(file, "# type: matrix\n");
  fprintf(file, "# rows: %d\n", rows);
  fprintf(file, "# columns: %d\n", cols);
}

}  // namespace octave
