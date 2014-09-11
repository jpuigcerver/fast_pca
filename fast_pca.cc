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
      "Usage: %s [-d] [-e eigval] [-f format] [-g eigvec] [-m mean]\n"
      "  [-o output] [-p idim] [-q odim] [-s stddev] [input ...]\n"
      "Options:\n"
      "  -d         use double precision\n"
      "  -e eigval  matrix with the eigenvalues of the covariance\n"
      "  -f format  format of the data matrix\n"
      "  -g eigvec  matrix with the eigenvectors of the covariance\n"
      "  -m mean    matrix with the mean of the data\n"
      "  -o output  output \n"
      "  -p idim    input dimensions\n"
      "  -q odim    output dimensions\n"
      "  -s stddev  matrix with the standard deviation of the data\n\n",
      prog);
}

template <typename real_t>
void load_pca_and_project(
    const string& format, const string& eigval_fn, const string& eigvec_fn,
    const string& mean_fn, const string& stdv_fn, int odim,
    const vector<string>& input, const string& output) {
  vector<real_t> eigval;
  vector<real_t> eigvec;
  vector<real_t> mean;
  vector<real_t> stdev;
  int one = 1, idim = -1;
  real_t cum_energy = -1;

  // load mean vector
  load_matlab<real_t>(mean_fn.c_str(), &one, &idim, &mean);
  // check input and output dimensions
  if (idim < odim) {
    fprintf(stderr, "ERROR: Invalid output dimension (%d > %d)!\n", odim, idim);
    exit(1);
  }
  // load eigenvectors
  load_matlab<real_t>(eigvec_fn.c_str(), &idim, &idim, &eigvec);
  // compute preserved cumulative energy
  if (eigval_fn != "") {
    load_matlab<real_t>(eigval_fn.c_str(), &one, &idim, &eigval);
    cum_energy = 0;
    for (int d = 0; d < odim; ++d) { cum_energy += eigval[d]; }
  }
  // load standard deviations
  if (stdv_fn != "") {
    load_matlab<real_t>(stdv_fn.c_str(), &one, &idim, &stdev);
  }
  // open input files
  vector<FILE*> input_files;
  if (input.size() == 0) { input_files.push_back(stdin); }
  for (size_t f = 0; f < input.size(); ++f) {
    FILE* file = fopen(input[f].c_str(), format == "binary" ? "rb" : "r");
    if (!file) {
      fprintf(stderr, "ERROR: Failed to open file \"%s\"!\n", input[f]);
      exit(1);
    }
    input_files.push_back(file);
  }

  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const string& fname = (file == stdin ? "**stdin**" : input[f]);
    bool ascii = true;
    if (format == "mat") {
      int num_rows = -1;
      read_matlab_header(fname, file, &num_rows, &idim);
    } else if (format == "binary") {
      ascii = false;
    }

    int tcols = 0;
    vector<real_t> x(idim);
    while (1) {
      const int tcols = ascii ?
          read_row<true, real_t>(file, idim, x) :
          read_row<false, real_t>(file, idim, x);
      if (tcols == 0) break;
      else if (tcols != idim) {
        fprintf(stderr, "ERROR: Corrupted matrix in file \"%s\"!\n", fname);
        exit(1);
      }
      project

    }

  }
}

int main(int argc, char** argv) {
  int opt = -1;
  int inp_dim = -1, out_dim = -1;
  bool simple = true;     // use simple precision ?
  string format = "mat";
  string eigval_fn = "";
  string eigvec_fn = "";
  string mean_fn = "";
  string stdv_fn = "";
  string output = "";
  while ((opt = getopt(argc, argv, "de:f:g:hm:n:o:s:")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'e':
        eigval_fn = optarg;
        break;
      case 'f':
        format = optarg;
        if (format != "mat" && format != "ascii" && format != "binary") {
          fprintf(stderr, "ERROR: Wrong matrix format: \"%s\"!\n", format.c_str());
          exit(1);
        }
        break;
      case 'g':
        eigvec_fn = optarg;
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'm':
        mean_fn = optarg;
        break;
      case 'n':
        inp_dim = atoi(optarg);
        if (inp_dim < 2) {
          fprintf(stderr, "ERROR: Input dimension must be greater than 1!\n");
          exit(1);
        }
        break;
      case 'o':
        output = optarg;
        break;
      case 's':
        stdv_fn = optarg;
        break;
      default:
        return 1;
    }
  }

  if (optind + 1 > argc) {
    fprintf(stderr, "ERROR: Missing arguments!\n\n");
    help(argv[0]);
    exit(1);
  }

  fprintf(stderr, "------------------------ Command line ------------------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (!simple) fprintf(stderr, " -d");
  if (eigval_fn != "") fprintf(stderr, " -e \"%s\"", eigval_fn.c_str());
  if (format != "") fprintf(stderr, " -f \"%s\"", format.c_str());
  if (eigvec_fn != "") fprintf(stderr, " -g \"%s\"", eigvec_fn.c_str());
  if (mean_fn != "") fprintf(stderr, " -m \"%s\"", mean_fn.c_str());
  if (inp_dim > 1) fprintf(stderr, " -n %d", inp_dim);
  if (output != "") fprintf(stderr, " -o \"%s\"", output.c_str());
  if (stdv_fn != "") fprintf(stderr, " -s \"%s\"", stdv_fn.c_str());
  for (int a = optind; a < argc - 1; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, " %s\n", argv[argc - 1]);
  fprintf(stderr, "--------------------------------------------------------------\n");

  // input file names
  vector<string> input;
  for(int a = optind; a < argc - 1; ++a) {
    input.push_back(argv[a]);
  }

  if (eigvec_fn != "" && mean_fn != "" && stdv_fn != "") {
    if (simple) {
      load_pca_and_project<float>(
          format, eigval_fn, eigvec_fn, mean_fn, stdv_fn, out_dim, output);
    } else {
      load_pca_and_project<double>(
          format, eigval_fn, eigvec_fn, mean_fn, stdv_fn, out_dim, output);
    }
    // load pca data and project
  } else if (input.size() > 0) {
    // compute pca data and project

  } else {
    fprintf(
        stderr,
        "ERROR: You must specify either -g, -m and -s options or an input file!\n");
    exit(1);
  }

}
