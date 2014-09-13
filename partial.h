#ifndef PARTIAL_H_
#define PARTIAL_H_

#include "math.h"
#include "file.h"

#include <numeric>
#include <string>
#include <vector>

using std::accumulate;
using std::string;
using std::vector;

// Compute sufficient variables to compute the covariance matrix
template <bool ascii, typename real_t>
int compute_partial(FILE* file, int dims, real_t* C, real_t* S) {
  vector<real_t> x(dims);
  int d = 0, n = 0;
  for (; (d = read_row<ascii, real_t>(file, dims, x.data())) == dims; ++n) {
    axpy<real_t>(dims, 1.0, x.data(), S);
    ger<real_t>(dims, dims, 1.0, x.data(), x.data(), C);
  }

  return (d == 0 || d == dims) ? n : -1;
}

template <typename real_t>
int reduce_partial(int mappers, int dims, int* n, real_t* C, real_t* S) {
  const int N = accumulate(n, n + mappers, 0);
  for (int m = 1; m < mappers; ++m) {
    axpy<real_t>(dims, 1.0, S + m * dims, S);
    axpy<real_t>(dims * dims, 1.0, C + m * dims * dims, C);
  }
  return N;
}

template <typename real_t>
int partial_cov_mean(
    const string& format, int *dims, const vector<string>& input, const string& output,
    vector<real_t>* m, vector<real_t>* c) {
  const bool ascii = format == "binary" ? false : true;
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
  vector<int> processed_rows(input_files.size(), 0);
  // read input headers
  if (format == "text") {
    for (size_t f = 0; f < input_files.size(); ++f) {
      FILE* file = input_files[f];
      const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
      read_text_header(fname, file, &expect_rows[f], dims);
    }
  } else if (*dims < 1) {
    fprintf(stderr, "ERROR: You must specify the input dimensions!\n");
    exit(1);
  }
  // this vector accumulates the sum across each dimension, needed to compute
  // the mean and the covariance of the data.
  m->resize(*dims, 0);
  // this matrix accumulates the sum of the cross-product across dimensions,
  // needed to compute the covariance of the data.
  c->resize((*dims) * (*dims), 0);
  // TODO(jpuigcerver): this can run in parallel
  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
    // compute partial results
    processed_rows[f] = ascii ?
        compute_partial<true, real_t>(
            file, *dims, c->data(), m->data()) :
        compute_partial<false, real_t>(
            file, *dims, c->data(), m->data());
    if (processed_rows[f] < 0 ||
        (expect_rows[f] > 0 && processed_rows[f] != expect_rows[f])) {
      fprintf(stderr, "ERROR: Corrupted matrix file \"%s\"!\n", fname);
      exit(1);
    }
  }
  return reduce_partial<real_t>(
      input_files.size(), *dims, processed_rows.data(), c->data(), m->data());
}

#endif  // PARTIAL_H_
