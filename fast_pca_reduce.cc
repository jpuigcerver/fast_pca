#include <getopt.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "partial.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s -d dims [-s] [input ...]\n"
      "Options:\n"
      "  -d dims     number of dimensions\n"
      "  -s          use simple precision (float numbers)\n",
      prog);
}

template <typename real_t>
void process_reduce_pca(const vector<string>& input, int dims) {
  vector<real_t> S(2 * dims, 0);
  vector<real_t> C(2 * dims * dims, 0);
  size_t N = 0;

  // accumulate partial results to compute the covariance matrix
  for (size_t f = 0; f < input.size(); ++f) {
    // load number of data samples
    int n = 0;
    const string num_fn = input[f] + ".num.part";
    load_integers<true>(num_fn.c_str(), 1, &n);
    if (n < 1) continue;
    // load partial covariance
    const string cov_fn = input[f] + ".cov.part";
    load_matrix<true, real_t>(
        cov_fn.c_str(), dims, dims, C.data() + dims * dims);
    // load partial mean
    const string mean_fn = input[f] + ".mean.part";
    load_matrix<true, real_t>(mean_fn.c_str(), 1, dims, S.data() + dims);
    // sum up partial results
    axpy<real_t>(dims, 1.0, S.data() + dims, 1, S.data(), 1);
    axpy<real_t>(dims * dims, 1.0, C.data() + dims * dims, 1, C.data(), 1);
    N += n;
  }

  // compute covariance matrix
  real_t* cm = C.data();
  for (int i = 0; i < dims; ++i) {
    for (int j = 0; j < dims; ++j) {
      cm[i * dims + j] = (cm[i * dims + j] - S[i] * S[j]) / (N - 1);
    }
  }

  dump_matrix<true, real_t>(stdout, dims, dims, cm);
}

int main(int argc, char** argv) {
  int opt = -1;
  int dims = 0;            // number of dimensions
  bool simple = false;     // use simple precision ?

  while ((opt = getopt(argc, argv, "d:sh")) != -1) {
    switch (opt) {
      case 'd':
        dims = atoi(optarg);
        break;
      case 's':
        simple = true;
        break;
      case 'h':
        help(argv[0]);
        return 0;
      default:
        return 1;
    }
  }

  fprintf(stderr, "-------------- Command line --------------\n");
  fprintf(stderr, "%s -d %d", argv[0], dims);
  if (simple) fprintf(stderr, " -s");
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "------------------------------------------\n");

  vector<string> input;
  for(int a = optind; a < argc; ++a) {
    input.push_back(argv[a]);
  }

  if (simple) {
    process_reduce_pca<float>(input, dims);
  } else {
    process_reduce_pca<double>(input, dims);
  }

  return 0;
}
