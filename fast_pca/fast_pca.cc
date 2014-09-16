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
#include "fast_pca/partial.h"
#include "fast_pca/pca.h"

using std::accumulate;
using std::min;
using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-C] [-P] [-d] [-s] [-p idim] [-q odim] [-j energy] "
      "[-f format] [-m pca] [-o output] [input ...]\n"
      "Options:\n"
      "  -C         compute pca from data\n"
      "  -P         project data using computed pca\n"
      "  -d         use double precision\n"
      "  -s         don't normalize data\n"
      "  -p idim    data input dimensions\n"
      "  -q odim    data output dimensions\n"
      "  -j energy  minimum cumulative energy preserved\n"
      "  -f format  format of the data matrix (ascii, binary, text)\n"
      "  -m pca     file containing the pca information\n"
      "  -o output  output data matrix\n",
      prog);
}

template <typename real_t>
void compute_pca(
    const string& format, const vector<string>& input, int* dims,
    vector<real_t>* eigval, vector<real_t>* eigvec, vector<real_t>* mean,
    vector<real_t>* stdev) {
  const bool ascii = (format == "binary" ? false : true);
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
  vector<int> process_rows(input_files.size(), 0);
  // read input headers
  if (format == "simple") {
    for (size_t f = 0; f < input_files.size(); ++f) {
      FILE* file = input_files[f];
      const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
      simple::read_header_ascii(fname, file, &expect_rows[f], dims);
    }
  } else if (*dims < 1) {
    fprintf(stderr, "ERROR: You must specify the input dimensions!\n");
    exit(1);
  }
  // Allocate memory for the results.
  // Take into account that the "mean" and "eigvec" will store partial results.
  mean->resize(*dims, 0);
  stdev->resize(*dims, 0);
  eigval->resize(*dims, 0);
  eigvec->resize(*dims * *dims, 0);
  // TODO(jpuigcerver): this can run in parallel
  for (size_t f = 0; f < input_files.size(); ++f) {
    FILE* file = input_files[f];
    const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
    // compute partial results
    process_rows[f] = ascii ?
        compute_partial<true, real_t>(
            file, (*dims), eigvec->data(), mean->data()) :
        compute_partial<false, real_t>(
            file, (*dims), eigvec->data(), mean->data());
    if (process_rows[f] < 0 ||
        (expect_rows[f] > 0 && process_rows[f] != expect_rows[f])) {
      fprintf(stderr, "ERROR: Corrupted matrix file \"%s\"!\n", fname);
      exit(1);
    }
  }
  // total number of processed rows
  const int N = accumulate(process_rows.begin(), process_rows.end(), 0);
  if (N < 2) {
    fprintf(stderr, "ERROR: Number of data rows is too small!\n");
    exit(1);
  }
  // compute pca
  if (pca<real_t>(
          N, *dims, mean->data(), eigvec->data(), stdev->data(),
          eigval->data()) != 0) {
    fprintf(stderr, "ERROR: Failed to compute PCA!\n");
    exit(1);
  }
  // compute energy quantiles
  // (how many dimensions you need to preserve % of energy)
  const double quant_val[] = {0.25, 0.5, 0.75, 1.0};
  int quant_dim[] = {*dims, *dims, *dims, *dims};
  const double total_energy = accumulate(eigval->begin(), eigval->end(), 0.0);
  double cum_energy = 0.0;
  for (int i = 0, q = 0; i < *dims && q < 4; ++i) {
    cum_energy += (*eigval)[i];
    for (; q < 4 && cum_energy >= total_energy * quant_val[q]; ++q) {
      quant_dim[q] = i + 1;
    }
  }
  fprintf(stderr, "-------------------- PCA summary --------------------\n");
  fprintf(stderr, "Processed rows: %d\n", N);
  fprintf(stderr, "Input dimension: %d\n", *dims);
  fprintf(
      stderr,
      "Energy quantiles: 25%% -> %d, 50%% -> %d, 75%% -> %d, 100%% -> %d\n",
      quant_dim[0], quant_dim[1], quant_dim[2], quant_dim[3]);
  fprintf(stderr, "-----------------------------------------------------\n");
}

template <typename real_t>
void project_data(
    const string& format,
    const vector<real_t>& eigval, const vector<real_t>& eigvec,
    const vector<real_t>& mean, const vector<real_t>& stdev,
    const vector<string>& input, const string& output, int idim, int odim,
    double min_energy, bool normalize_data) {
  const bool ascii = (format == "binary" ? false : true);
  // compute cumulative energy preserved and (optionally) output dimension
  const double total_energy = accumulate(eigval.begin(), eigval.end(), 0.0);
  double cum_energy = -1.0;
  if (odim > 0) {
    cum_energy =
        accumulate(eigval.begin(), eigval.begin() + odim, 0.0) / total_energy;
  } else {
    cum_energy = eigval[0];
    for (odim = 1; odim < idim && (cum_energy / total_energy < min_energy);
         ++odim) {
      cum_energy += eigval[odim];
    }
    cum_energy = min(cum_energy / total_energy, 1.0);
  }
  // check input and output dimensions
  if (idim < odim) {
    fprintf(stderr, "ERROR: Invalid output dimension (%d > %d)!\n", odim, idim);
    exit(1);
  }
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
  if (format == "simple") {
    // input headers
    int total_rows = 0;
    for (size_t f = 0; f < input_files.size(); ++f) {
      FILE* file = input_files[f];
      const char* fname = (file == stdin ? "**stdin**" : input[f].c_str());
      int file_rows = -1;
      simple::read_header_ascii(fname, file, &file_rows, &idim);
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
      if (project(
              1, idim, odim, eigvec.data(), mean.data(),
              normalize_data ? stdev.data() : NULL, x.data()) != 0) {
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
  fprintf(stderr, "---------------- Projection summary -----------------\n");
  fprintf(stderr, "Processed rows: %d\n", processed_rows);
  fprintf(stderr, "Input dimension: %d\n", idim);
  fprintf(stderr, "Output dimension: %d\n", odim);
  fprintf(stderr, "Preserved energy: %.4g%%\n", cum_energy * 100.0);
  fprintf(stderr, "-----------------------------------------------------\n");
}

template <typename real_t>
void do_work(
    const bool do_compute_pca, const bool do_project_data,
    const string& format, const string& pca_fn, const vector<string>& input,
    const string& output, int inp_dim, int out_dim, double min_energy,
    bool normalize_data) {
  vector<real_t> mean;
  vector<real_t> stdev;
  vector<real_t> eigval;
  vector<real_t> eigvec;
  if (do_compute_pca) {
    compute_pca<real_t>(
        format, input, &inp_dim, &eigval, &eigvec, &mean, &stdev);
    if (pca_fn != "") {
      save_pca<real_t>(pca_fn.c_str(), inp_dim, mean, stdev, eigval, eigvec);
    } else if (!do_project_data || output != "") {
      save_pca<real_t>(NULL, inp_dim, mean, stdev, eigval, eigvec);
    }
  } else {
    if (pca_fn == "") {
      fprintf(stderr, "ERROR: Projection requires pca file (-m)!\n");
      exit(1);
    }
    load_pca<real_t>(pca_fn.c_str(), &inp_dim, &mean, &stdev, &eigval, &eigvec);
  }
  if (do_project_data) {
    project_data<real_t>(
        format, eigval, eigvec, mean, stdev, input, output, inp_dim, out_dim,
        min_energy, normalize_data);
  }
}

int main(int argc, char** argv) {
  int opt = -1;
  int inp_dim = -1, out_dim = -1;
  bool simple_precision = true;     // use simple precision ?
  bool normalize_data = true;
  string format = "simple";
  string pca_fn = "";
  string output = "";
  double min_energy = -1.0;
  bool do_compute_pca = false;
  bool do_project_data = false;
  while ((opt = getopt(argc, argv, "CPdshf:j:m:o:p:q:")) != -1) {
    switch (opt) {
      case 'C':
        do_compute_pca = true;
        break;
      case 'P':
        do_project_data = true;
        break;
      case 'd':
        simple_precision = false;
        break;
      case 's':
        normalize_data = false;
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'f':
        format = optarg;
        if (format != "simple" && format != "ascii" && format != "binary") {
          fprintf(
              stderr, "ERROR: Wrong matrix format: \"%s\"!\n", format.c_str());
          exit(1);
        }
        break;
      case 'j':
        min_energy = atof(optarg);
        if (min_energy <= 0.0 || min_energy > 1.0) {
          fprintf(
              stderr, "ERROR: Invalid minimum cumulative energy: %f!\n",
              min_energy);
        }
        break;
      case 'm':
        pca_fn = optarg;
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
      default:
        return 1;
    }
  }

  if (!do_compute_pca && !do_project_data) {
    do_compute_pca = do_project_data = true;
  }

  fprintf(stderr, "-------------------- Command line -------------------\n");
  fprintf(stderr, "%s", argv[0]);
  if (do_compute_pca) fprintf(stderr, " -C");
  if (do_project_data) fprintf(stderr, " -P");
  if (!simple_precision) fprintf(stderr, " -d");
  if (!normalize_data) fprintf(stderr, " -s");
  if (inp_dim > 0) fprintf(stderr, " -p %d", inp_dim);
  if (out_dim > 0) fprintf(stderr, " -q %d", out_dim);
  if (min_energy > 0) fprintf(stderr, " -j %g", min_energy);
  if (format != "") fprintf(stderr, " -f \"%s\"", format.c_str());
  if (pca_fn != "") fprintf(stderr, " -m \"%s\"", pca_fn.c_str());
  if (output != "") fprintf(stderr, " -o \"%s\"", output.c_str());
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  if (min_energy < 0) min_energy = 1.0;

  // input file names
  vector<string> input;
  for (int a = optind; a < argc; ++a) {
    input.push_back(argv[a]);
  }
  // when reading from stdin, we cannot process data twice!
  if (do_compute_pca && do_project_data && input.size() == 0) {
    fprintf(stderr, "ERROR: When reading from stdin, use either -C or -P!\n");
    exit(1);
  }
  if (simple_precision) {
    do_work<float>(
        do_compute_pca, do_project_data, format, pca_fn, input, output,
        inp_dim, out_dim, min_energy, normalize_data);
  } else {
    do_work<double>(
        do_compute_pca, do_project_data, format, pca_fn, input, output,
        inp_dim, out_dim, min_energy, normalize_data);
  }

  return 0;
}
