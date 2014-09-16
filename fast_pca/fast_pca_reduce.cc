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
#include "fast_pca/pca.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-d] input [input ...] output\n"
      "Options:\n"
      "  -d         use double precision\n",
      prog);
}

template <typename real_t>
void do_work(const vector<string>& input, const string& output) {
  vector<real_t> m;
  vector<real_t> c;
  vector<real_t> mean;
  vector<real_t> eigvec;
  int processed_rows = 0;
  int dims = -1;  // number of dimensions will be read from MAT file

  // accumulate partial results to compute the covariance matrix
  for (size_t f = 0; f < input.size(); ++f) {
    // load partial results
    int n = -1;
    partial::load_partial<real_t>(input[f].c_str(), &n, &dims, &m, &c);
    // sum up partial results
    if (mean.size() == 0) { mean.resize(dims); }
    if (eigvec.size() == 0) { eigvec.resize(dims * dims); }
    axpy<real_t>(dims, 1.0, m.data(), mean.data());
    axpy<real_t>(dims * dims, 1.0, c.data(), eigvec.data());
    processed_rows += n;
  }
  vector<real_t> stdev(dims);
  vector<real_t> eigval(dims);
  if (pca<real_t>(
          processed_rows, dims, mean.data(), eigvec.data(), stdev.data(),
          eigval.data()) != 0) {
    fprintf(stderr, "ERROR: Failed to compute pca!\n");
    exit(1);
  }
  save_pca<real_t>(output.c_str(), dims, mean, stdev, eigval, eigvec);
}

int main(int argc, char** argv) {
  int opt = -1;
  bool simple = true;     // use simple precision ?
  while ((opt = getopt(argc, argv, "de:hs:")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
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
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  if (optind + 2 > argc) {
    fprintf(stderr, "ERROR: Missing arguments!\n");
    exit(1);
  }

  vector<string> input;
  for (int a = optind; a < argc - 1; ++a) {
    input.push_back(argv[a]);
  }

  if (simple) {
    do_work<float>(input, argv[argc - 1]);
  } else {
    do_work<double>(input, argv[argc - 1]);
  }

  return 0;
}
