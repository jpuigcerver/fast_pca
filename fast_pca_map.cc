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
      "Usage: %s [-d] [-p dim] [-f format] [input ...] out_prefix\n"
      "Options:\n"
      "  -d          use double precision\n"
      "  -p dim      data dimensions\n"
      "  -f format   format of the data matrix (ascii, binary, text)\n",
      prog);
}

template <typename real_t>
void do_work(
    const string& format, int dims, const vector<string>& input,
    const string& output) {
  vector<real_t> m;
  vector<real_t> c;
  const int n = partial_cov_mean<real_t>(format, &dims, input, output, &m, &c);
  save_integers<true>((output + ".rows.part").c_str(), 1, &n);
  save_text<real_t>((output + ".mean.part").c_str(), 1, dims, m.data());
  save_text<real_t>((output + ".cov.part").c_str(), dims, dims, c.data());
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

  if (simple) {
    do_work<float>(format, dims, input, output);
  } else {
    do_work<double>(format, dims, input, output);
  }

  return 0;
}
