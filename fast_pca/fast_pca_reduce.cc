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
#include "fast_pca/file_pca.h"
#include "fast_pca/pca.h"


using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-d] [-o output] [input ...]\n"
      "Options:\n"
      "  -d         use double precision\n"
      "  -o output  output file\n",
      prog);
}

template <typename real_t>
void do_work(const vector<string>& input, const string& output) {
  vector<real_t> M;  // global mean
  vector<real_t> C;  // global co-momentum
  vector<real_t> D;  // diff between the global mean and the file mean
  vector<real_t> m;
  vector<real_t> c;
  int n = -1;
  int dims = -1;
  // process first file
  load_n_mean_cov<real_t>(input[0], &n, &dims, &M, &C);
  D.resize(dims);
  // process rest of files
  for (size_t f = 1; f < input.size(); ++f) {
    // load n, dims, mean and covariance
    int br = -1;
    load_n_mean_cov<real_t>(input[f], &br, &dims, &m, &c);
    // D = M - m
    memcpy(D.data(), M.data(), sizeof(real_t) * dims);
    axpy<real_t>(dims, -1, m.data(), D.data());
    // update co-moments matrix
    // C += c
    axpy<real_t>(dims * dims, 1, c.data(), C.data());
    // C += D * D' * (br * n) / (br + n)
    ger<real_t>(
        dims, dims, (1.0 * br) * n / (n + br), D.data(), D.data(), C.data());
    // update mean
    for (int d = 0; d < dims; ++d) {
      M[d] = (n * M[d] + br * m[d]) / (n + br);
    }
    n += br;
  }
  LOG("%d", n);
  // convert comoment into covariance matrix
  for (int i = 0; i < dims * dims; ++i) { C[i] /= (n - 1); }
  // compute standard deviation in each dimension
  vector<real_t> stddev(dims);
  for (int i = 0; i < dims; ++i) { stddev[i] = sqrt(C[i * dims + i]); }
  // compute eigenvectors and eigenvalues of the covariance matrix
  // WARNING: This destroys the covariance matrix!
  vector<real_t> eigval(dims);
  eig<real_t>(dims, C.data(), eigval.data());
  save_pca(output, dims, M, stddev, eigval, C);
  /*
  vector<real_t> stdev(dims);
  vector<real_t> eigval(dims);
  if (pca<real_t>(
          processed_rows, dims, mean.data(), eigvec.data(), stdev.data(),
          eigval.data()) != 0) {
    fprintf(stderr, "ERROR: Failed to compute pca!\n");
    exit(1);
  }
  save_pca<real_t>(output.c_str(), dims, mean, stdev, eigval, eigvec);*/
}

int main(int argc, char** argv) {
  int opt = -1;
  bool simple = true;     // use simple precision ?
  string output = "";     // output filename
  while ((opt = getopt(argc, argv, "dho:")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'o':
        output = optarg;
        break;
      default:
        return 1;
    }
  }

  fprintf(stderr, "-------------------- Command line -------------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (!simple) fprintf(stderr, " -d");
  if (output != "") fprintf(stderr, " -o \"%s\"", output.c_str());
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  vector<string> input;
  for (int a = optind; a < argc; ++a) { input.push_back(argv[a]); }
  if (input.empty()) input.push_back("");

  if (simple) {
    do_work<float>(input, output);
  } else {
    do_work<double>(input, output);
  }

  return 0;
}
