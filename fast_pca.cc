#include <getopt.h>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "file.h"
#include "pca.h"
#include "partial.h"

using std::accumulate;
using std::min;
using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-d] [-e eigval] [-f format] [-g eigvec] [-j energy] "
      "[-m mean] [-o output] [-p idim] [-q odim] [-s stddev] [input ...]\n"
      "Options:\n"
      "  -d         use double precision\n"
      "  -e eigval  matrix with the eigenvalues of the covariance\n"
      "  -f format  format of the data matrix\n"
      "  -g eigvec  matrix with the eigenvectors of the covariance\n"
      "  -j energy  minimum cumulative energy preserved\n"
      "  -m mean    matrix with the mean of the data\n"
      "  -o output  output \n"
      "  -p idim    input dimensions\n"
      "  -q odim    output dimensions\n"
      "  -s stddev  matrix with the standard deviation of the data\n\n",
      prog);
}

template <typename real_t>
void compute_pca(
    const string& format,
    const vector<string>& input, int* dims,
    vector<real_t>* eigval, vector<real_t>* eigvec,
    vector<real_t>* mean, vector<real_t>* stdev) {
  const bool ascii = (format == "binary" ? false : true);
  // reserve memory
  vector<real_t> C((*dims) * (*dims), 0);
  vector<real_t> S((*dims), 0);
  vector<int> n(input.size(), 0);
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
  // TODO(jpuigcerver): this can run in parallel
  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
    int expected_rows = -1;
    if (format == "mat") {
      read_matlab_header(fname, file, &expected_rows, dims);
    }
    n[f] = ascii ?
        compute_partial<true, real_t>(
            file, (*dims), C.data() + f * (*dims) * (*dims),
            S.data() + f * (*dims)) :
        compute_partial<false, real_t>(
            file, (*dims), C.data() + f * (*dims) * (*dims),
            S.data() + f * (*dims));
    if (n[f] < 0 ||
        (expected_rows > 0 && n[f] != expected_rows)) {
      fprintf(stderr, "ERROR: Corrupted matrix file \"%s\"!\n", fname);
      exit(1);
    }
  }
  // reduce data
  const int N = reduce_partial<real_t>(
      input.size(), (*dims), n.data(), C.data(), S.data());

  // compute means
  for (int i = 0; i < (*dims); ++i) {
    mean[i] /= N;
  }
  // compute covariance matrix
  for (int i = 0; i < (*dims); ++i) {
    for (int j = 0; j < (*dims); ++j) {
      cov[i * (*dims) + j] =
          (cov[i * (*dims) + j] - N * mean[i] * mean[j]) / (N - 1);
    }
  }
  // compute eigenvalues of the covariance matrix (cov matrix is destroyed)
  real_t* w = part_m.data();
  vector<real_t> s((*dims));
  if (pca<real_t>((*dims), cov.data(), w, s.data()) != 0) {
    fprintf(stderr, "ERROR: Failed to compute PCA!\n");
    exit(1);
  }
}

template <typename real_t>
double load_pca(
    const string& eigval_fn, const string& eigvec_fn, const string& mean_fn,
    const string& stdv_fn, int *odim, double min_energy,
    vector<real_t>* eigval, vector<real_t>* eigvec,
    vector<real_t>* mean, vector<real_t>* stdev) {
  double cum_energy = -1;
  int idim = -1, one = 1;
  // load mean vector
  load_matlab<real_t>(mean_fn.c_str(), &one, &idim, mean);
  // check input and output dimensions
  if (idim < *odim) {
    fprintf(
        stderr, "ERROR: Invalid output dimension (%d > %d)!\n", *odim, idim);
    exit(1);
  }
  // load eigenvectors
  load_matlab<real_t>(eigvec_fn.c_str(), &idim, &idim, eigvec);
  // compute preserved cumulative energy
  if (eigval_fn != "") {
    load_matlab<real_t>(eigval_fn.c_str(), &one, &idim, eigval);
    const double total_energy = accumulate(eigval->begin(), eigval->end(), 0);
    if (*odim > 0)  {
      // if the output dimension is given, compute the preserved energy
      cum_energy = accumulate(eigval->begin(), eigval->begin() + *odim, 0) /
          total_energy;
    } else {
      // if the output dimension is not given, ensure at least min_energy is
      // preserved.
      cum_energy = (*eigval)[0];
      for (*odim = 1; *odim < idim && (cum_energy / total_energy < min_energy);
           ++(*odim)) {
        cum_energy += (*eigval)[*odim];
      }
      cum_energy = min(cum_energy / total_energy, 1.0);
    }
  } else if (*odim < 1) {
    *odim = idim;
  }
  // load standard deviations
  if (stdv_fn != "") {
    load_matlab<real_t>(stdv_fn.c_str(), &one, &idim, stdev);
  }
  return cum_energy;
}

template <typename real_t>
int project_data(
    const string& format,
    const vector<real_t>& eigval, const vector<real_t>& eigvec,
    const vector<real_t>& mean, const vector<real_t>& stdev,
    const int odim, const vector<string>& input,
    const string& output) {
  const bool ascii = (format == "binary" ? false : true);
  int idim = mean.size();
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
  // open output file
  FILE* output_file = stdout;
  if (output != "" && !(output_file = fopen(
          output.c_str(), format == "binary" ? "wb" : "w"))) {
    fprintf(stderr, "ERROR: Failed to open file \"%s\"!\n", output.c_str());
    exit(1);
  }
  // process MAT headers
  if (format == "mat") {
    // input headers
    int total_rows = 0;
    for (size_t f = 0; f < input_files.size(); ++f) {
      FILE* file = input_files[f];
      const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
      int file_rows = -1;
      read_matlab_header(fname, file, &file_rows, &idim);
      total_rows += file_rows;
    }
    // output header
    fprintf(output_file, "%d %d\n", total_rows, odim);
  }
  // process input files
  int processed_rows = 0;
  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
    int tcols = 0;
    vector<real_t> x(idim);
    while (1) {
      // read data row
      const int tcols = ascii ?
          read_row<true, real_t>(file, idim, x.data()) :
          read_row<false, real_t>(file, idim, x.data());
      if (tcols == 0) break;
      else if (tcols != idim) {
        fprintf(stderr, "ERROR: Corrupted matrix in file \"%s\"!\n", fname);
        exit(1);
      }
      // project data row
      if (project_single(
              idim, odim, eigvec.data(), mean.data(),
              stdev.size() ? stdev.data() : NULL, x.data()) != 0) {
        fprintf(
            stderr, "ERROR: Projection failed in file \"%s\"!\n", fname);
        exit(1);
      }
      // output projected data
      if (ascii) {
        write_row<true, real_t>(output_file, odim, x.data());
      } else {
        write_row<false, real_t>(output_file, odim, x.data());
      }
      ++processed_rows;
    }
    fclose(file);
  }
  return processed_rows;
}


template <typename real_t>
void load_pca_and_project(
    const string& format, const string& eigval_fn, const string& eigvec_fn,
    const string& mean_fn, const string& stdv_fn, int odim, double min_energy,
    const vector<string>& input, const string& output) {
  vector<real_t> eigval;
  vector<real_t> eigvec;
  vector<real_t> mean;
  vector<real_t> stdev;
  const double cum_energy = load_pca<real_t>(
      eigval_fn, eigvec_fn, mean_fn, stdv_fn, &odim, min_energy,
      &eigval, &eigvec, &mean, &stdev);
  const int processed_rows = project_data<real_t>(
      format, eigval, eigvec, mean, stdev, odim, input, output);
  fprintf(stderr, "---------------------- Summary ----------------------\n");
  fprintf(stderr, "Processed rows: %d\n", processed_rows);
  fprintf(stderr, "Input dimension: %d\n", static_cast<int>(mean.size()));
  fprintf(stderr, "Output dimension: %d\n", odim);
  fprintf(stderr, "Cumulative energy: %g\n", cum_energy);
  fprintf(stderr, "-----------------------------------------------------\n");
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
  double min_energy = 1.0;
  while ((opt = getopt(argc, argv, "de:f:g:j:hm:o:p:q:s:")) != -1) {
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
          fprintf(
              stderr, "ERROR: Wrong matrix format: \"%s\"!\n", format.c_str());
          exit(1);
        }
        break;
      case 'g':
        eigvec_fn = optarg;
        break;
      case 'j':
        min_energy = atof(optarg);
        if (min_energy <= 0.0 || min_energy > 1.0) {
          fprintf(
              stderr, "ERROR: Invalid minimum cumulative energy: %f!\n",
              min_energy);
        }
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'm':
        mean_fn = optarg;
        break;
      case 'o':
        output = optarg;
        break;
      case 'p':
        inp_dim = atoi(optarg);
        if (inp_dim < 1) {
          fprintf(stderr, "ERROR: Input dimension must be positive!\n");
          exit(1);
        }
        break;
      case 'q':
        out_dim = atoi(optarg);
        if (out_dim < 1) {
          fprintf(stderr, "ERROR: Output dimension must be positive!\n");
          exit(1);
        }
        break;
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
  if (eigval_fn != "") fprintf(stderr, " -e \"%s\"", eigval_fn.c_str());
  if (format != "") fprintf(stderr, " -f \"%s\"", format.c_str());
  if (eigvec_fn != "") fprintf(stderr, " -g \"%s\"", eigvec_fn.c_str());
  if (mean_fn != "") fprintf(stderr, " -m \"%s\"", mean_fn.c_str());
  if (output != "") fprintf(stderr, " -o \"%s\"", output.c_str());
  if (inp_dim > 0) fprintf(stderr, " -p %d", inp_dim);
  if (out_dim > 0) fprintf(stderr, " -q %d", out_dim);
  if (stdv_fn != "") fprintf(stderr, " -s \"%s\"", stdv_fn.c_str());
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  // input file names
  vector<string> input;
  for(int a = optind; a < argc; ++a) {
    input.push_back(argv[a]);
  }

  if (eigvec_fn != "" && mean_fn != "") {
    // load pca data and project
    if (simple) {
      load_pca_and_project<float>(
          format, eigval_fn, eigvec_fn, mean_fn, stdv_fn, out_dim, min_energy,
          input, output);
    } else {
      load_pca_and_project<double>(
          format, eigval_fn, eigvec_fn, mean_fn, stdv_fn, out_dim, min_energy,
          input, output);
    }
  } else if (input.size() > 0) {
    // compute pca data and project
    if (simple) {
      //compute_pca_and_project<float>(
      //    format, inp_dim, out_dim, min_energy, input, output);
    } else {
    }

  } else {
    fprintf(
        stderr,
        "ERROR: You must specify either -g, -m and -s options or an input!\n");
    exit(1);
  }

}
