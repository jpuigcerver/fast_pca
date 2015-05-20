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
using std::max;
using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-C] [-P] [-b size] [-d] [-e dim] [-f format] [-j energy] "
      "[-m pca] [-n] [-o output] [-p idim] [-q odim] [input ...]\n"
      "Options:\n"
      "  -C         compute pca from data\n"
      "  -P         project data using computed pca\n"
      "  -b size    process data in batches of this number of rows\n"
      "  -d         use double precision\n"
      "  -e dim     exclude first (positive) or last (negative) dimensions\n"
      "  -f format  format of the data matrix (ascii, binary, octave, vbosch)\n"
      "  -j energy  minimum cumulative energy preserved\n"
      "  -m pca     file containing the pca information\n"
      "  -n         normalize data before projection\n"
      "  -o output  output data matrix\n"
      "  -p idim    data input dimensions\n"
      "  -q odim    data output dimensions\n",
      prog);
}


template <FORMAT_CODE fmt, typename real_t>
void compute_pca(
    vector<string> input, int block, int exclude_dims, int* dims,
    vector<real_t>* eigval, vector<real_t>* eigvec, vector<real_t>* mean,
    vector<real_t>* stddev) {
  int n = 0;  // number of data samples
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
  // Up to here, we have the whole covariance matrix. But we will compute
  // only the eigenvalues and eigenvectors of the submatrix obtained after
  // removing the first/last `exclude_dims' rows and columns.
  const int eff_dims = *dims - abs(exclude_dims);
  real_t* eigvec_offset = exclude_dims <= 0 ? eigvec->data() :
      eigvec->data() + exclude_dims * (*dims) + exclude_dims;
  // compute eigenvectors and eigenvalues of the covariance matrix
  // WARNING: This destroys the covariance matrix!
  eigval->resize(eff_dims);
  CHECK(eig<real_t>(eff_dims, *dims, eigvec_offset, eigval->data()) == 0);
  // Move all eigenvectors to the first rows
  for (int r = 0; r < eff_dims; ++r) {
    for (int d = 0; d < eff_dims; ++d) {
      (*eigvec)[r * eff_dims + d] = eigvec_offset[r * (*dims) + d];
    }
  }
  eigvec->resize(eff_dims * eff_dims);
}

template <FORMAT_CODE fmt, typename real_t>
void project_data(
    vector<string> input, const string& output, int block, int idim, int odim,
    int exclude_dims, double min_energy, bool normalize_data,
    const vector<real_t>& mean, const vector<real_t>& stddev,
    const vector<real_t>& eigval, const vector<real_t>& eigvec) {
  CHECK_MSG(idim > 0, "Invalid input dimension!");
  CHECK_MSG(odim > 0 || min_energy > 0.0,
            "Specify output dimensions or minimum energy!");
  CHECK_MSG(odim < 0 || odim >= exclude_dims,
        "Number of output dimensions (%d) is lower than the number of "
        "excluded dimensions from PCA (%d)!", odim, exclude_dims);
  // compute cumulative energy preserved and (optionally) output dimension
  const double total_energy = accumulate(eigval.begin(), eigval.end(), 0.0);
  double cum_energy = 0.0;
  if (odim > 0) {
    cum_energy = accumulate(
        eigval.begin(), eigval.begin() + odim - abs(exclude_dims), 0.0);
    cum_energy /= total_energy;
  } else {
    cum_energy = 0.0;
    for (odim = 0; odim < idim - abs(exclude_dims) &&
             (cum_energy / total_energy < min_energy); ++odim) {
      cum_energy += eigval[odim];
    }
    odim += abs(exclude_dims);
    odim = max(odim, 1);
    cum_energy = min(cum_energy / total_energy, 1.0);
  }
  // check input and output dimensions
  CHECK_MSG(odim <= idim, "Invalid output dimension!");
  // open input files
  vector<FILE*> files;
  open_files("rb", "**stdin**", stdin, &input, &files);
  // process input headers
  int expected_rows = -1;
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = input[f].c_str();
    FILE* file = files[f];
    int h_rows = -1, h_dims = -1;
    CHECK_MSG(
        !read_matrix_header<fmt>(file, NULL, &h_rows, &h_dims),
        "Invalid header in file \"%s\"!", fname);
    CHECK_MSG(
        (h_dims < 1 || idim == h_dims),
        "Bad number of dimensions in file \"%s\"!", fname);
    if (expected_rows < 0) expected_rows = h_rows;
    else expected_rows += h_rows;
  }
  FILE* output_file = stdout;
  if (output != "") { output_file = open_file(output.c_str(), "wb"); }
  write_matrix_header<fmt>(output_file, "", expected_rows, odim);
  // ----- process input files -----
  vector<real_t> x(block * idim, 0);  // data block
  vector<real_t> b(block * idim, 0);  // auxiliar data block
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
          br, idim, odim, exclude_dims, eigvec.data(), mean.data(),
          normalize_data ? stddev.data() : NULL, x.data(), b.data());
      // output the computed data
      write_matrix<fmt, real_t>(output_file, br, odim, x.data());
      // update total number of processed rows
      n += br;
    }
  }
  if (output != "") { fclose(output_file); }
  CHECK_MSG(
      (expected_rows < 1 || expected_rows == n),
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
void pca_summary(int inp_dims, int exclude_dims, const vector<real_t>& eigval) {
  const int eff_dims = inp_dims - abs(exclude_dims);
  const double total_energy = accumulate(eigval.begin(), eigval.end(), 0.0);
  // compute energy quantiles
  // (how many dimensions you need to preserve % of energy)
  const double quant_val[] = {0.25, 0.5, 0.75, 1.0};
  int quant_dim[] = {eff_dims, eff_dims, eff_dims, eff_dims};
  double cum_energy = 0.0;
  for (int i = 0, q = 0; i < eff_dims && q < 4; ++i) {
    cum_energy += eigval[i];
    for (; q < 4 && cum_energy >= total_energy * quant_val[q]; ++q) {
      quant_dim[q] = i + 1;
    }
  }
  fprintf(
      stderr,
      "-------------------- PCA summary --------------------\n"
      "Input dimensions: %d\n"
      "Excluded dimensions: %d\n"
      "Energy quantiles: 25%% -> %d, 50%% -> %d, 75%% -> %d, 100%% -> %d\n"
      "-----------------------------------------------------\n",
      inp_dims, exclude_dims, quant_dim[0], quant_dim[1], quant_dim[2],
      quant_dim[3]);
}

template <FORMAT_CODE fmt, typename real_t>
void do_work(
    const bool do_compute_pca, const bool do_project_data,
    const string& pca_fn, const vector<string>& input,
    const string& output, int block, int inp_dim, int out_dim,
    double min_energy, bool normalize_data, int exclude_dims) {
  vector<real_t> mean;
  vector<real_t> stdev;
  vector<real_t> eigval;
  vector<real_t> eigvec;
  if (do_compute_pca) {
    compute_pca<fmt, real_t>(
        input, block, exclude_dims, &inp_dim, &eigval, &eigvec, &mean, &stdev);
    if (!do_project_data || pca_fn != "" || output != "") {
      save_pca<real_t>(
          pca_fn, inp_dim, exclude_dims, mean, stdev, eigval, eigvec);
    }
  } else {
    CHECK_MSG(pca_fn != "", "Specify a pca file to load from!");
    if (exclude_dims != 0) {
      WARN("Ignoring \"-e %d\": excluded dimensions will be read from the "
           "pca file...", exclude_dims);
    }
    load_pca<real_t>(
        pca_fn, &inp_dim, &exclude_dims, &mean, &stdev, &eigval, &eigvec);
  }
  pca_summary(inp_dim, exclude_dims, eigval);
  if (do_project_data) {
    project_data<fmt, real_t>(
        input, output, block, inp_dim, out_dim, exclude_dims, min_energy,
        normalize_data, mean, stdev, eigval, eigvec);
  }
}

int main(int argc, char** argv) {
  int opt = -1;
  int inp_dim = -1, out_dim = -1;
  int exclude_dims = 0;
  int block = 1000;
  bool simple_precision = true;     // use simple precision ?
  bool normalize_data = false;
  string pca_fn = "";
  string output = "";
  double min_energy = -1.0;
  bool do_compute_pca = false;
  bool do_project_data = false;
  FORMAT_CODE format = FMT_ASCII;
  const char* format_str = NULL;
  while ((opt = getopt(argc, argv, "CPb:de:f:hj:m:no:p:q:")) != -1) {
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
      case 'e':
        exclude_dims = atoi(optarg);
        break;
      case 'n':
        normalize_data = true;
        break;
      case 'f':
        format_str = optarg;
        format = format_code_from_name(format_str);
        CHECK_MSG(format != FMT_UNKNOWN, "Unknown format (-f \"%s\")!", optarg);
        break;
      case 'h':
        help(argv[0]);
        return 0;
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
  if (exclude_dims) fprintf(stderr, " -e %d", exclude_dims);
  if (normalize_data) fprintf(stderr, " -n");
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

  // Launch the appropiate do_work function, depending on the format of the
  // data and whether double or single precision is used.
  // NOTE: The reason for the extreme `if-else' branching here, is due to the
  // fact that `do_work' is a templated function. Actually, during compile-time
  // several `do_work' instances are compiled, and here we must call the correct
  // one.
  switch (format) {
    case FMT_ASCII:
      if (simple_precision) {
        do_work<FMT_ASCII, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      } else {
        do_work<FMT_ASCII, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      }
      break;
    case FMT_BINARY:
      if (simple_precision) {
        do_work<FMT_OCTAVE, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      } else {
        do_work<FMT_OCTAVE, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      }
      break;
    case FMT_OCTAVE:
      if (simple_precision) {
        do_work<FMT_OCTAVE, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      } else {
        do_work<FMT_OCTAVE, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      }
      break;
    case FMT_VBOSCH:
      if (simple_precision) {
        do_work<FMT_VBOSCH, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      } else {
        do_work<FMT_VBOSCH, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      }
      break;
    default:
      ERROR("Not implemented for this format!");
  }

  return 0;
}
