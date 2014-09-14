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
      "Usage: %s [-d] [-e eigval] [-s stddev] input [input ...] means eigvec\n"
      "Options:\n"
      "  -d         use double precision\n"
      "  -e eigval  save eigenvalues to this file\n"
      "  -s stddev  save the std. deviations to this file\n",
      prog);
}

template <typename real_t>
void do_work(
    const vector<string>& input, const string& mean_fn, const string& eigvec_fn,
    const string& eigval_fn, const string& stdev_fn) {
  vector<real_t> m;
  vector<real_t> c;
  vector<real_t> mean;
  vector<real_t> eigvec;
  int processed_rows = 0;
  int dims = -1;  // number of dimensions will be read from MAT file
  int one = 1;    // used to read the mean partial files

  // accumulate partial results to compute the covariance matrix
  for (size_t f = 0; f < input.size(); ++f) {
    // load number of data samples
    int n = -1;
    load_integers<true>((input[f] + ".rows.part").c_str(), 1, &n);
    if (n < 1) continue;
    // load partial mean
    load_text<real_t>((input[f] + ".mean.part").c_str(), &one, &dims, &m);
    // load partial covariance
    load_text<real_t>((input[f] + ".cov.part").c_str(), &dims, &dims, &c);
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
  save_text<real_t>(mean_fn.c_str(), 1, dims, mean.data());
  save_text<real_t>(eigvec_fn.c_str(), dims, dims, eigvec.data());
  if (eigval_fn != "") {
    save_text<real_t>(eigval_fn.c_str(), 1, dims, eigval.data());
  }
  if (stdev_fn != "") {
    save_text<real_t>(stdev_fn.c_str(), 1, dims, stdev.data());
  }
}

int main(int argc, char** argv) {
  int opt = -1;
  bool simple = true;     // use simple precision ?
  string eigv_fn = "";
  string stdv_fn = "";
  while ((opt = getopt(argc, argv, "de:hs:")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'e':
        eigv_fn = optarg;
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 's':
        stdv_fn = optarg;
        break;
      default:
        return 1;
    }
  }

  fprintf(stderr, "-------------------- Command line -------------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (!simple) fprintf(stderr, " -d");
  if (eigv_fn != "") fprintf(stderr, " -e \"%s\"", eigv_fn.c_str());
  if (stdv_fn != "") fprintf(stderr, " -s \"%s\"", stdv_fn.c_str());
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  if (optind + 3 > argc) {
    fprintf(stderr, "ERROR: Missing input arguments!\n");
    exit(1);
  }

  vector<string> input;
  for (int a = optind; a < argc - 2; ++a) {
    input.push_back(argv[a]);
  }

  if (simple) {
    do_work<float>(input, argv[argc - 2], argv[argc - 1], eigv_fn, stdv_fn);
  } else {
    do_work<double>(input, argv[argc - 2], argv[argc - 1], eigv_fn, stdv_fn);
  }

  return 0;
}
