/*
  The MIT License (MIT)

  Copyright (c) 2014 Joan Puigcerver

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <getopt.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/partial.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-d] [-p dim] [-f format] [-o output] [input ...]\n"
      "Options:\n"
      "  -d          use double precision\n"
      "  -p dim      data dimensions\n"
      "  -f format   format of the data matrix (ascii, binary, text)\n"
      "  -o output   output file\n",
      prog);
}

template <typename real_t>
void do_work(
    const string& format, int block_rows, int dims, const string& output,
    vector<string>* input) {
  if (format == "plain" && dims <= 0) {
    fprintf(stderr, "ERROR: You must specify the number of dimensions!\n");
    exit(1);
  }
  // open input files
  vector<FILE*> files;
  open_files(
      format == "binary" ? "rb" : "r", "**stdin**", stdin, input, &files);
  // read header information
  vector<int> expect_rows(files->size(), -1);
  for (size_t f = 0; f < files.size(); ++f) {
    if (format == "simple") {
      simple::read_header_ascii(
          input[f].c_str(), files[f], &expect_rows[f], &dims);
    } else if (format == "octave") {
      // TODO
      exit(2);
    }
  }
  vector<real_t> x(block_rows * dims);
  vector<real_t> m(dims);
  vector<real_t> c(dims * dims);
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    read_block<real_t>(file, block_rows * dims, x->data());
  }


}

int main(int argc, char** argv) {
  int opt = -1;
  int dims = -1;             // number of dimensions
  bool simple = true;        // use simple precision ?
  string format = "simple";  // input matrix format
  string output = "";

  while ((opt = getopt(argc, argv, "df:o:p:h")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'f':
        format = optarg;
        if (format != "simple" && format != "ascii" && format != "binary") {
          fprintf(stderr, "ERROR: Unknown format!\n");
          exit(1);
        }
        break;
      case 'o':
        output = optarg;
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

  fprintf(stderr, "-------------------- Command line -------------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (!simple) fprintf(stderr, " -d");
  if (dims > 0) fprintf(stderr, " -p %d", dims);
  fprintf(stderr, " -f \"%s\"", format.c_str());
  if (output != "") fprintf(stderr, "-o %s", output.c_str());
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  vector<string> input;
  for (int a = optind; a < argc; ++a) { input.push_back(argv[a]); }

  if (simple) {
    do_work<float>(format, dims, input, output);
  } else {
    do_work<double>(format, dims, input, output);
  }

  return 0;
}
