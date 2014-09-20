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

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/pca.h"
#include "fast_pca/fast_pca.h"
#include "fast_pca/logging.h"

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
      "  -b size    process data in batches of this number of rows\n"
      "  -d         use double precision\n"
      "  -f format  format of the data matrix (ascii, binary, text)\n"
      "  -j energy  minimum cumulative energy preserved\n"
      "  -m pca     file containing the pca information\n"
      "  -n         don't normalize data\n"
      "  -o output  output data matrix\n"
      "  -p idim    data input dimensions\n"
      "  -q odim    data output dimensions\n",
      prog);
}


template <FORMAT_CODE fmt, typename real_t>
void compute_pca(
    vector<string> input, int block, int* dims,
    vector<real_t>* eigval, vector<real_t>* eigvec, vector<real_t>* mean,
    vector<real_t>* stddev) {
  int n = 0;
  // process input to compute mean and co-moments
  compute_mean_comoments_from_inputs<fmt, real_t>(
      block, input, &n, dims, mean, eigvec);
  // compute covariance from co-moments
  CHECK_MSG(n > 1, "You need at least 2 data points (%d processed)!", n);
  for (int i = 0; i < (*dims) * (*dims); ++i) {
    (*eigvec)[i] /= (n - 1);
  }
  // compute standard deviation in each dimension
  stddev->resize(*dims);
  for (int i = 0; i < (*dims); ++i) {
    (*stddev)[i] = sqrt((*eigvec)[i * (*dims) + i]);
  }
  // compute eigenvectors and eigenvalues of the covariance matrix
  // WARNING: This destroys the covariance matrix!
  eigval->resize(*dims);
  CHECK(eig<real_t>(*dims, eigvec->data(), eigval->data()) == 0);
}

template <FORMAT_CODE fmt, typename real_t>
void project_data(
    vector<string> input, const string& output, int block, int idim, int odim,
    int min_energy, bool normalize_data, const vector<real_t>& mean,
    const vector<real_t>& stddev, const vector<real_t>& eigval,
    const vector<real_t>& eigvec) {
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
  CHECK_MSG(odim <= idim, "Invalid output dimension!");
  // open input files
  vector<FILE*> files;
  open_files("rb", "**stdin**", stdin, &input, &files);
  // process input headers
  int expected_rows = 0;
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int h_rows = -1, h_dims = -1;
    CHECK_MSG(
        !read_matrix_header<fmt>(file, NULL, &h_rows, &h_dims),
        "Invalid header in file \"%s\"!", fname);
    CHECK_MSG(
        idim == h_dims, "Bad number of dimensions in file \"%s\"!", fname);
    expected_rows += h_rows;
  }
  FILE* output_file = stdout;
  if (output != "") { output_file = open_file(output.c_str(), "wb"); }
  write_matrix_header<fmt>(output_file, "", expected_rows, odim);
  // ----- process input files -----
  vector<real_t> x(block * idim, 0);  // data block
  int n = 0;  // total processed rows
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int fr = 0, be = 0, br = 0;
    while ((be = read_block<fmt, real_t>(file, block * idim, x.data())) > 0) {
      CHECK_MSG(be % idim == 0, "Corrupted matrix in file \"%s\"!", fname);
      br = be / idim;
      fr += br;
      // project input data using pca. WARNING: destroys the original data!
      project<real_t>(
          br, idim, odim, eigvec.data(), mean.data(),
          normalize_data ? stddev.data() : NULL, x.data());
      // output the computed data
      write_matrix<fmt, real_t>(output_file, br, odim, x.data());
      // update total number of processed rows
      n += br;
    }
  }
  if (output != "") { fclose(output_file); }
  CHECK_MSG(
      expected_rows == n,
      "Number of processed rows is lower than expected "
      "(expected: %d, processed: %d)!", expected_rows, n);
  fprintf(stderr, "---------------- Projection summary -----------------\n");
  fprintf(stderr, "Processed rows: %d\n", n);
  fprintf(stderr, "Input dimension: %d\n", idim);
  fprintf(stderr, "Output dimension: %d\n", odim);
  fprintf(stderr, "Preserved energy: %.4g%%\n", cum_energy * 100.0);
  fprintf(stderr, "-----------------------------------------------------\n");
}

template <typename real_t>
void pca_summary(int dims, const vector<real_t>& eigval) {
  const double total_energy = accumulate(eigval.begin(), eigval.end(), 0.0);
  // compute energy quantiles
  // (how many dimensions you need to preserve % of energy)
  const double quant_val[] = {0.25, 0.5, 0.75, 1.0};
  int quant_dim[] = {dims, dims, dims, dims};
  double cum_energy = 0.0;
  for (int i = 0, q = 0; i < dims && q < 4; ++i) {
    cum_energy += eigval[i];
    for (; q < 4 && cum_energy >= total_energy * quant_val[q]; ++q) {
      quant_dim[q] = i + 1;
    }
  }
  fprintf(
      stderr,
      "-------------------- PCA summary --------------------\n"
      "Input dimension: %d\n"
      "Energy quantiles: 25%% -> %d, 50%% -> %d, 75%% -> %d, 100%% -> %d\n"
      "-----------------------------------------------------\n",
      dims, quant_dim[0], quant_dim[1], quant_dim[2], quant_dim[3]);
}

template <FORMAT_CODE fmt, typename real_t>
void do_work(
    const bool do_compute_pca, const bool do_project_data,
    const string& pca_fn, const vector<string>& input,
    const string& output, int block, int inp_dim, int out_dim,
    double min_energy, bool normalize_data) {
  vector<real_t> mean;
  vector<real_t> stdev;
  vector<real_t> eigval;
  vector<real_t> eigvec;
  if (do_compute_pca) {
    compute_pca<fmt, real_t>(
        input, block, &inp_dim, &eigval, &eigvec, &mean, &stdev);
    if (!do_project_data || pca_fn != "" || output != "") {
      save_pca<real_t>(pca_fn, inp_dim, mean, stdev, eigval, eigvec);
    }
  } else {
    CHECK_MSG(pca_fn != "", "Specify a pca file to load from!");
    load_pca<real_t>(pca_fn, &inp_dim, &mean, &stdev, &eigval, &eigvec);
  }
  pca_summary(inp_dim, eigval);
  if (do_project_data) {
    project_data<fmt, real_t>(
        input, output, block, inp_dim, out_dim, min_energy, normalize_data,
        mean, stdev, eigval, eigvec);
  }
}

int main(int argc, char** argv) {
  int opt = -1;
  int inp_dim = -1, out_dim = -1;
  int block = 1000;
  bool simple_precision = true;     // use simple precision ?
  bool normalize_data = true;
  string pca_fn = "";
  string output = "";
  double min_energy = -1.0;
  bool do_compute_pca = false;
  bool do_project_data = false;
  FORMAT_CODE format = FMT_VBOSCH;
  const char* format_str = NULL;
  while ((opt = getopt(argc, argv, "CPb:df:hj:m:no:p:q:")) != -1) {
    switch (opt) {
      case 'C':
        do_compute_pca = true;
        break;
      case 'P':
        do_project_data = true;
        break;
      case 'b':
        block = atoi(optarg);
        CHECK_MSG(block > 0, "Block size must be positive (-b %d)!", block);
        break;
      case 'd':
        simple_precision = false;
        break;
      case 'n':
        normalize_data = false;
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'f':
        format_str = optarg;
        format = format_code_from_name(format_str);
        CHECK_MSG(format != FMT_UNKNOWN, "Unknown format (-f \"%s\")!", optarg);
        break;
      case 'j':
        min_energy = atof(optarg);
        CHECK_MSG(
            min_energy >= 0.0 && min_energy <= 1.0,
            "Invalid minimum cumulative energy (-j %f)!", min_energy);
        break;
      case 'm':
        pca_fn = optarg;
        break;
      case 'o':
        output = optarg;
        break;
      case 'p':
        inp_dim = atoi(optarg);
        CHECK_MSG(
            inp_dim > 0, "Input dimension must be positive (-p %d)!", inp_dim);
        break;
      case 'q':
        out_dim = atoi(optarg);
        CHECK_MSG(
            out_dim > 0, "Output dimension must be positive (-q %d)!", out_dim);
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
  if (format_str) fprintf(stderr, " -f \"%s\"", format_str);
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
  CHECK_MSG(
      !do_compute_pca || !do_project_data || input.size() > 0,
      "Use either -C or -P when reading from stdin!");

  switch (format) {
    case FMT_ASCII:
      if (simple_precision) {
        do_work<FMT_ASCII, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      } else {
        do_work<FMT_ASCII, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      }
      break;
    case FMT_BINARY:
      if (simple_precision) {
        do_work<FMT_OCTAVE, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      } else {
        do_work<FMT_OCTAVE, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      }
      break;
    case FMT_OCTAVE:
      if (simple_precision) {
        do_work<FMT_OCTAVE, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      } else {
        do_work<FMT_OCTAVE, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      }
      break;
    case FMT_VBOSCH:
      if (simple_precision) {
        do_work<FMT_VBOSCH, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      } else {
        do_work<FMT_VBOSCH, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data);
      }
      break;
    default:
      ERROR("Not implemented for this format!");
  }

  return 0;
}
