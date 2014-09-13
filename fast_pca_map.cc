#include <getopt.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "file.h"
#include "partial.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-d] [-p dim] [-f format] [input ...] output\n"
      "Options:\n"
      "  -d          use double precision\n"
      "  -p dim      data dimensions\n"
      "  -f format   format of the data matrix (ascii, binary, text)\n",
      prog);
}

int main(int argc, char** argv) {
  int opt = -1;
  int dims = -1;          // number of dimensions
  string format = "text"; // input matrix format
  bool simple = true;     // use simple precision ?

  while ((opt = getopt(argc, argv, "d:f:o:sh")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'f':
        format = optarg;
        break;
      case 'p':
        dims = atoi(optarg);
        if (dims < 1) {
          fprintf(stderr, "ERROR: Input dimension must be positive!\n");
          exit(1);
        }
        break;
      case 'h':
        help(argv[0]);
        return 0;
      default:
        return 1;
    }
  }

  fprintf(stderr, "-------------- Command line --------------\n");
  if (!simple) fprintf(stderr, " -d");
  if (format != "") fprintf(stderr, " -f \"%s\"", format.c_str());
  if (dims > 0) fprintf(stderr, " -p %d", dims);
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "------------------------------------------\n");

  if (format != "text" && format != "ascii" && format != "binary") {
    fprintf(stderr, "ERROR: Unknown format!\n");
    exit(1);
  }

  if (optind + 1 > argc) {
    fprintf(stderr, "ERROR: You must specify an output prefix!\n");
    exit(1);
  }

  const string& output = argv[argc - 1];
  vector<string> input;
  for (int a = optind; a < argc - 1; ++a) { input.push_back(argv[a]); }

  // open input files
  vector<FILE*> input_files;
  if (input.size() == 0) { input_files.push_back(stdin); }
  for (size_t f = 0; f < input.size(); ++f) {
    FILE* file = fopen(input[f].c_str(), format == "binary" ? "rb" : "r");
    if (!file) {
      fprintf(stderr, "ERROR: Failed to open file \"%s\"!\n", input[f].c_str());
      exit(1);
    }
    input_files.push_back(file);
  }
  // number of expected rows in each input file
  vector<int> expect_rows(input_files.size(), -1);
  // number of processed rows in each input file
  vector<int> process_rows(input_files.size(), 0);
  // read input headers
  if (format == "text") {
    for (size_t f = 0; f < input_files.size(); ++f) {
      FILE* file = input_files[f];
      const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
      read_text_header(fname, file, &expect_rows[f], &dims);
    }
  } else if (dims < 1) {
    fprintf(stderr, "ERROR: You must specify the input dimensions!\n");
    exit(1);
  }
  // Allocate memory for the results.
  //
  // ALLOCATE DATA
  //
  // TODO(jpuigcerver): this can run in parallel
  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
    // compute partial results
    process_rows[f] = ascii ?
        compute_partial<true, real_t>(
            file, (*dims), eigvec->data(), mean->data()) :
        compute_partial<false, real_t>(
            file, (*dims), eigvec->data(), mean->data());
    if (process_rows[f] < 0 ||
        (expect_rows[f] > 0 && process_rows[f] != expect_rows[f])) {
      fprintf(stderr, "ERROR: Corrupted matrix file \"%s\"!\n", fname);
      exit(1);
    }
  }



  return 0;
}
