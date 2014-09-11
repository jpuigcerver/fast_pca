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
      "Usage: %s [-d dims] [-f format] [-s] -o output [input ...]\n"
      "Options:\n"
      "  -d dims     number of dimensions\n"
      "  -f format   format of the input matrix (ascii, binary, mat)\n"
      "  -o output   prefix of the output files\n"
      "  -s          use simple precision (float numbers)\n\n"
      "Matrix formats:\n"
      " - ascii)\n"
      "  A matrix is represented by a string of float/double numbers\n"
      "  represented using your current locale.\n"
      " - binary)\n"
      "  A matrix is represented by a chunk of data containing float/double\n"
      "  numbers. This is dependent of the byte order of your machine!\n"
      " - mat)\n"
      "  A matrix is represented using the ASCII MAT file format.\n"
      "  The first line contains the number of rows and columns of the\n"
      "  matrix, and each new line contains a row with the elements separated\n"
      "  by whitespaces.\n\n"
      "Important considerations:\n"
      " - The default format is \"ascii\".\n"
      " - All matrices are assumed to be in Row-major order (C-order).\n"
      " - The number of dimensions is optional only for \"mat\" format.\n",
      prog);
}

// -----------------------------------------------------
// Process MAT files
// -----------------------------------------------------
template <typename real_t>
void process_mat_files(
    char** argv, const int optind, const vector<FILE*>& input,
    const string& output) {
  int dims = 0;       // number of dimensions
  vector<int> n;      // number of data points
  n.push_back(0);
  // read matrix dimensions from the first input file
  if (fscanf(input[0], "%d %d", &n.back(), &dims) != 2) {
    fprintf(
        stderr, "ERROR: Failed to read header in \"%s\"!\n", argv[optind]);
    exit(1);
  }
  // check whether the matrix has a valid size
  if (dims == 0 || n.back() == 0) {
    fprintf(stderr, "ERROR: Zero-size matrix in \"%s\"!\n", argv[optind]);
    exit(1);
  }
  // read header from other files and check it...
  for (size_t f = 1; f < input.size(); ++f) {
    int d = 0;
    n.push_back(0);
    // read header
    if (fscanf(input[f], "%d %d", &n.back(), &d) != 2) {
      fprintf(
          stderr,
          "ERROR: Failed to read header in \"%s\"!n", argv[optind + f]);
      exit(1);
    }
    // check whether the dimensions are consistent
    if (d != dims) {
      fprintf(
          stderr, "ERROR: Invalid dimension size in \"%s\"!\n",
          argv[optind + f]);
      exit(1);
    }
    // check for a valid data size
    if (n.back() == 0) {
      fprintf(
          stderr, "ERROR: Zero-size matrix in \"%s\"!\n", argv[optind + f]);
      exit(1);
    }
  }

  // reserve memory
  vector<real_t> C(dims * dims, 0);
  vector<real_t> S(dims, 0);

  // this can run in parallel
  for (size_t f = 0; f < input.size(); ++f) {
    const char* fname = (input[f] == stdin ? "*stdin*" : argv[optind + f]);
    fprintf(stderr, "File: \"%s\", Rows: %d, Cols: %d\n", fname, n[f], dims);
    if (compute_partial<true, real_t>(
            input[f], dims, C.data() + f * dims * dims,
            S.data() + f * dims) != n[f]) {
      fprintf(stderr, "ERROR: Bad MAT file \"%s\"!\n", fname);
      exit(1);
    }
  }

  // reduce data
  const int N = reduce_partial<real_t>(
      input.size(), dims, n.data(), C.data(), S.data());

  // number of processed samples
  const string num_fn = output + ".num.part";
  save_integers<true>(num_fn.c_str(), 1, &N);
  // partial covariance
  const string cov_fn = output + ".cov.part";
  save_matrix<true, real_t>(cov_fn.c_str(), dims, dims, C.data());
  // partial mean
  const string mean_fn = output + ".mean.part";
  save_matrix<true, real_t>(mean_fn.c_str(), 1, dims, S.data());
}


// -----------------------------------------------------
// Process simple files
// -----------------------------------------------------
template <bool ascii, typename real_t>
void process_simple_files(
    char** argv, const int optind, const int dims, const vector<FILE*>& input,
    const string& output) {
  vector<real_t> C(dims * dims);
  vector<real_t> S(dims);
  vector<int> n(input.size());

  // this can run in parallel
  for (size_t f = 0; f < input.size(); ++f) {
    n[f] = compute_partial<ascii, real_t>(
        input[f], dims, C.data() + f * dims * dims,
        S.data() + f * dims) != n[f];
    if (n[f] < 0) {
      fprintf(stderr, "ERROR: Bad matrix file \"%s\"!\n", argv[optind + f]);
      exit(1);
    }
  }

  // reduce data
  const int N = reduce_partial<real_t>(
      input.size(), dims, n.data(), C.data(), S.data());

  // number of processed samples
  const string num_fn = output + ".num.part";
  save_integers<true>(num_fn.c_str(), 1, &N);
  // partial covariance
  const string cov_fn = output + ".cov.part";
  save_matrix<true, real_t>(cov_fn.c_str(), dims, dims, C.data());
  // partial mean
  const string mean_fn = output + ".mean.part";
  save_matrix<true, real_t>(mean_fn.c_str(), 1, dims, S.data());
}

int main(int argc, char** argv) {
  int opt = -1;
  int dims = 0;            // number of dimensions
  string format = "ascii"; // input matrix format
  string output = "";      // output path prefix
  bool simple = false;     // use simple precision ?

  while ((opt = getopt(argc, argv, "d:f:o:sh")) != -1) {
    switch (opt) {
      case 'd':
        dims = atoi(optarg);
        break;
      case 'f':
        format = optarg;
        break;
      case 'o':
        output = optarg;
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
  fprintf(stderr, "%s -d %d -f \"%s\" -o \"%s\"", argv[0], dims, format.c_str(),
          output.c_str());
  if (simple) fprintf(stderr, " -s");
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n");
  fprintf(stderr, "------------------------------------------\n");


  if (format != "ascii" && format != "binary" && format != "mat") {
    fprintf(stderr, "ERROR: Invalid format \"%s\"!\n", format.c_str());
    exit(1);
  }
  const bool mat_fmt = (format == "mat");
  const bool binary = (format == "binary");

  vector<FILE*> input;
  for(int a = optind; a < argc; ++a) {
    FILE* f = fopen(argv[a], binary ? "rb" : "r");
    if (!f) {
      fprintf(stderr, "ERROR: Failed to open file \"%s\"!\n", argv[a]);
      exit(1);
    }
    input.push_back(f);
  }
  if (input.size() == 0) {
    input.push_back(stdin);
  }

  if (!mat_fmt && dims < 1) {
    fprintf(stderr, "ERROR: Number of dimensions (-d) must be at least 1!\n");
    exit(1);
  }

  if (output == "") {
    fprintf(stderr, "ERROR: You must specify an output prefix (-o)!\n");
    exit(1);
  }

  if (mat_fmt) {
    if (simple) {
      process_mat_files<float>(argv, optind, input, output);
    } else {
      process_mat_files<double>(argv, optind, input, output);
    }
  } else if (binary) {
    if (simple) {
      process_simple_files<false, float>(argv, optind, dims, input, output);
    } else {
      process_simple_files<false, double>(argv, optind, dims, input, output);
    }
  } else {
    if (simple) {
      process_simple_files<true, float>(argv, optind, dims, input, output);
    } else {
      process_simple_files<true, double>(argv, optind, dims, input, output);
    }
  }

  return 0;
}
