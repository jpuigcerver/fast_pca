/*
  The MIT License (MIT)

  Copyright (c) 2014,2015 Joan Puigcerver

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

#include "fast_pca/fast_pca_common.h"
#include "fast_pca/file.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/pca.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [options] [input ...]\n\n"
      "Options:\n"
      "  -c         do not compute eigenvalues; output co-moments instead\n"
      "  -d         use double precision\n"
      "  -e dims    do not project first (positive) or last (negative) dims\n"
      "  -j energy  minimum relative amount of energy preserved\n"
      "  -m output  write (temporal) pca information to this file\n"
      "  -q odim    maximum output dimensions of the projected data\n",
      prog);
}

template <typename real_t>
void do_work(
    const vector<string>& input, const string& output, bool compute_pca,
    int exclude_dims, int out_dim, double min_rel_energy) {
  vector<real_t> M;  // global mean
  vector<real_t> C;  // global co-momentum
  vector<real_t> D;  // diff between the global mean and the file mean
  vector<real_t> m;
  vector<real_t> c;
  int n = -1;
  int inp_dim = -1;
  double miss_energy = 0.0;
  // process first file
  load_n_mean_cov<real_t>(input[0], &n, &inp_dim, &M, &C);
  D.resize(inp_dim);
  // process rest of files
  for (size_t f = 1; f < input.size(); ++f) {
    // load n, inp_dim, mean and covariance
    int br = -1;
    load_n_mean_cov<real_t>(input[f], &br, &inp_dim, &m, &c);
    // D = M - m
    memcpy(D.data(), M.data(), sizeof(real_t) * inp_dim);
    axpy<real_t>(inp_dim, -1, m.data(), D.data());
    // update co-moments matrix
    // C += c
    axpy<real_t>(inp_dim * inp_dim, 1, c.data(), C.data());
    // C += D * D' * (br * n) / (br + n)
    ger<real_t>(
        inp_dim, inp_dim, (1.0 * br) * n / (n + br),
        D.data(), D.data(), C.data());
    // update mean
    for (int d = 0; d < inp_dim; ++d) {
      M[d] = (n * M[d] + br * m[d]) / (n + br);
    }
    n += br;
  }
  if (compute_pca) {
    CHECK_FMT(
        inp_dim >= exclude_dims,
        "Dimensions to exclude (%d) is bigger than the data "
        "dimensions (%d)!", exclude_dims, inp_dim);
    CHECK_FMT(
        inp_dim >= out_dim,
        "Number of output dimensions (%d) is greater than the data "
        "dimensions (%d)!", out_dim, inp_dim);
    CHECK_FMT(
        out_dim < 1 || out_dim >= abs(exclude_dims),
        "Number of non-projected dimensions (%d) is bigger than the output "
        "dimensions (%d)!", abs(exclude_dims), out_dim);
    // convert comoment into covariance matrix
    for (int i = 0; i < inp_dim * inp_dim; ++i) { C[i] /= (n - 1); }
    // compute standard deviation in each dimension
    vector<real_t> stddev(inp_dim);
    for (int i = 0; i < inp_dim; ++i) { C[i] = sqrt(C[i * inp_dim + i]); }
    // compute eigenvectors and eigenvalues of the covariance matrix
    vector<real_t> eigval;
    compute_pca_from_covariance<real_t>(
        exclude_dims, min_rel_energy, inp_dim, &out_dim, &miss_energy,
        &C, &eigval);
    // Compute PCA summary
    vector<real_t> cumulative_energy;
    compute_cumulative_energy(eigval, &cumulative_energy);
    pca_summary<real_t>(
        inp_dim, exclude_dims, miss_energy, cumulative_energy);
    save_pca<real_t>(
        output, exclude_dims, miss_energy, M, stddev, eigval, C);
  } else {
    save_n_mean_cov(output, n, inp_dim, M, C);
  }
}

int main(int argc, char** argv) {
  int opt = -1;
  bool simple = true;        // use simple precision ?
  bool compute_pca = true;   // compute pca from the co-moments matrices
  int exclude_dims = 0;      // exclude these first/last dimensions from PCA
  string output = "";        // output filename
  int out_dim = -1;          // output dimension
  double min_rel_energy = -1.0;  // preserve energy
  while ((opt = getopt(argc, argv, "cde:hj:m:q:")) != -1) {
    switch (opt) {
      case 'c':
        compute_pca = false;
        break;
      case 'd':
        simple = false;
        break;
      case 'e':
        exclude_dims = atoi(optarg);
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'j':
        min_rel_energy = atof(optarg);
        CHECK_FMT(
            min_rel_energy >= 0.0 && min_rel_energy <= 1.0,
            "Invalid minimum amount of relative energy (-j %f)!",
            min_rel_energy);
        break;
      case 'm':
        output = optarg;
        break;
      case 'q':
        out_dim = atoi(optarg);
        CHECK_FMT(
            out_dim > 0, "Output dimension must be positive (-q %d)!", out_dim);
        break;
      default:
        return 1;
    }
  }

  fprintf(stderr, "-------------------- Command line -------------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (!compute_pca) fprintf(stderr, " -c");
  if (!simple) fprintf(stderr, " -d");
  if (exclude_dims) fprintf(stderr, " -e %d", exclude_dims);
  if (min_rel_energy > 0) fprintf(stderr, " -j %g", min_rel_energy);
  if (output != "") fprintf(stderr, " -m \"%s\"", output.c_str());
  if (out_dim > 0) fprintf(stderr, " -q %d", out_dim);
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  vector<string> input;
  for (int a = optind; a < argc; ++a) { input.push_back(argv[a]); }
  if (input.empty()) input.push_back("");

  if (simple) {
    do_work<float>(
        input, output, compute_pca, exclude_dims, out_dim, min_rel_energy);
  } else {
    do_work<double>(
        input, output, compute_pca, exclude_dims, out_dim, min_rel_energy);
  }

  return 0;
}
