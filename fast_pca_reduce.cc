#include <getopt.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "file.h"
#include "pca.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-d] [-e eigval] [-s stddev] input [input ...] means eigvec\n"
      "Options:\n"
      "  -d         use double precision\n"
      "  -e eigval  save eigenvalues to this file\n"
      "  -h         show help\n\n"
      "  -s stddev  save the std. deviations to this file\n"
      "This tool reads the partial results produced by fast_pca_partial\n"
      "to compute the covariance matrix of a dataset and computes the\n"
      "the mean vector and the eigenvectors of the covariance. Optionally,\n"
      "it can compute the eigenvalues of the covariance matrix as well.\n\n"
      "For each input matrix X, the program tries to load these files:\n"
      "   X.num.part  -> number of rows in X\n"
      "   X.cov.part  -> X * X'\n"
      "   X.mean.part -> Ones(1, D) * X\n",
      prog);
}

template <typename real_t>
void process_reduce_pca(
    const vector<string>& input, const string& mean_fn, const string& proj_fn,
    const string& eigv_fn, const string& stdv_fn) {
  vector<real_t> part_m;
  vector<real_t> part_c;
  vector<real_t> mean;
  vector<real_t> cov;
  size_t N = 0;
  int dims = -1;  // number of dimensions will be read from MAT file
  int one = 1;    // used to read the mean partial files

  // accumulate partial results to compute the covariance matrix
  for (size_t f = 0; f < input.size(); ++f) {
    // load number of data samples
    int n = 0;
    const string num_fn = input[f] + ".num.part";
    load_integers<true>(num_fn.c_str(), 1, &n);
    if (n < 1) continue;
    // load partial covariance
    const string cov_fn = input[f] + ".cov.part";
    load_matlab<real_t>(cov_fn.c_str(), &dims, &dims, &part_c);
    // load partial mean
    const string mean_fn = input[f] + ".mean.part";
    load_matlab<real_t>(mean_fn.c_str(), &one, &dims, &part_m);
    // sum up partial results
    if (mean.size() == 0) { mean.resize(dims); }
    if (cov.size() == 0) { cov.resize(dims * dims); }
    axpy<real_t>(dims, 1.0, part_m.data(), mean.data());
    axpy<real_t>(dims * dims, 1.0, part_c.data(), cov.data());
    N += n;
  }
  // compute means
  for (int i = 0; i < dims; ++i) {
    mean[i] /= N;
  }
  // compute covariance matrix
  for (int i = 0; i < dims; ++i) {
    for (int j = 0; j < dims; ++j) {
      cov[i * dims + j] =
          (cov[i * dims + j] - N * mean[i] * mean[j]) / (N - 1);
    }
  }
  // compute eigenvalues of the covariance matrix (cov matrix is destroyed)
  real_t* w = part_m.data();
  vector<real_t> s(dims);
  if (pca<real_t>(dims, cov.data(), w, s.data()) != 0) {
    fprintf(stderr, "ERROR: Failed to compute PCA!\n");
    exit(1);
  }


  save_matlab<real_t>(mean_fn.c_str(), 1, dims, mean.data());
  save_matlab<real_t>(proj_fn.c_str(), dims, dims, cov.data());
  if (eigv_fn != "") {
    save_matlab<real_t>(eigv_fn.c_str(), 1, dims, w);
  }
  if (stdv_fn != "") {
    save_matlab<real_t>(stdv_fn.c_str(), 1, dims, s.data());
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

  if (optind + 3 > argc) {
    fprintf(stderr, "ERROR: Missing arguments!\n\n");
    help(argv[0]);
    exit(1);
  }

  fprintf(stderr, "-------------- Command line --------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (!simple) fprintf(stderr, " -d");
  if (eigv_fn != "") fprintf(stderr, " -e \"%s\"", eigv_fn.c_str());
  if (stdv_fn != "") fprintf(stderr, " -s \"%s\"", stdv_fn.c_str());
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "------------------------------------------\n");

  vector<string> input;
  for(int a = optind; a < argc - 2; ++a) {
    input.push_back(argv[a]);
  }

  if (simple) {
    process_reduce_pca<float>(
        input, argv[argc - 2], argv[argc - 1], eigv_fn, stdv_fn);
  } else {
    process_reduce_pca<double>(
        input, argv[argc - 2], argv[argc - 1], eigv_fn, stdv_fn);
  }

  return 0;
}
