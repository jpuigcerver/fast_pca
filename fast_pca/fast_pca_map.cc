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
#include <cstring>
#include <string>
#include <vector>

#include "fast_pca/fast_pca.h"
#include "fast_pca/logging.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-b size] [-d] [-f format] [-o output] [-d dim] [input ...]\n"
      "Options:\n"
      "  -b size    process data in batches of this number of rows\n"
      "  -d         use double precision\n"
      "  -f format  format of the data matrix\n"
      "  -o output  output file\n"
      "  -p dim     data dimensions\n",
      prog);
}

template <FORMAT_CODE fmt, typename real_t>
void do_work(int block, int dims, string output, vector<string> input) {
  int n;
  vector<real_t> M;  // global mean
  vector<real_t> C;  // global co-moments matrix
  // compute mean and comoments matrix
  compute_mean_comoments_from_inputs<fmt, real_t>(
      block, input, &n, &dims, &M, &C);
  // output number of processed rows, mean and co-moments matrix
  save_n_mean_cov(output, n, dims, M, C);
}

int main(int argc, char** argv) {
  int opt = -1;
  int dims = -1;             // number of dimensions
  int block = 1000;          // block size
  bool simple = true;        // use simple precision ?
  string output = "";
  FORMAT_CODE format = FMT_ASCII;
  const char* format_str = NULL;

  while ((opt = getopt(argc, argv, "db:f:o:p:h")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'b':
        block = atoi(optarg);
        CHECK_MSG(block > 0, "Block size must be positive (-b %d)!", block);
        break;
      case 'f':
        format_str = optarg;
        format = format_code_from_name(format_str);
        CHECK_MSG(format != FMT_UNKNOWN, "Unknown format (-f \"%s\")!", optarg);
        break;
      case 'o':
        output = optarg;
        break;
      case 'p':
        dims = atoi(optarg);
        CHECK_MSG(dims > 0, "Input dimensions must be positive (-p %d)!", dims);
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
  fprintf(stderr, " -b %d", block);
  if (!simple) fprintf(stderr, " -d");
  if (format_str) fprintf(stderr, " -f \"%s\"", format_str);
  if (output != "") fprintf(stderr, "-o %s", output.c_str());
  if (dims > 0) fprintf(stderr, " -p %d", dims);
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  vector<string> input;
  for (int a = optind; a < argc; ++a) { input.push_back(argv[a]); }

  switch (format) {
    case FMT_ASCII:
      if (simple)
        do_work<FMT_ASCII, float>(block, dims, output, input);
      else
        do_work<FMT_ASCII, double>(block, dims, output, input);
      break;
    case FMT_BINARY:
      if (simple)
        do_work<FMT_BINARY, float>(block, dims, output, input);
      else
        do_work<FMT_BINARY, double>(block, dims, output, input);
      break;
    case FMT_OCTAVE:
      if (simple)
        do_work<FMT_OCTAVE, float>(block, dims, output, input);
      else
        do_work<FMT_OCTAVE, double>(block, dims, output, input);
      break;
    case FMT_VBOSCH:
      if (simple)
        do_work<FMT_VBOSCH, float>(block, dims, output, input);
      else
        do_work<FMT_VBOSCH, double>(block, dims, output, input);
      break;
    default:
      ERROR("Not implemented for this format!");
  }

  return 0;
}
