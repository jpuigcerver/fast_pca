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

#include <algorithm>
#include <numeric>
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
      "Usage: %s [-C] [-P] [options] [...]\n\n"
      "Examples:\n"
      "Compute PCA: %s -C [options] [input ...]\n"
      "Project: %s -P -m pca.mat [options] [input [output] ...]\n"
      "Compute PCA & project: %s -C -P [options] input [input ...] [output]\n\n"
      "Options:\n"
      "  -C         compute pca from data\n"
      "  -P         project data using computed pca\n"
      "  -b size    number of rows in the batch (default: 1000)\n"
      "  -d         use double precision\n"
      "  -e dim     exclude first (positive) or last (negative) dimensions\n"
      "  -f format  format of the data matrix (ascii, binary, octave, vbosch,\n"
      "             prhlt_htk)\n"
      "  -j energy  minimum cumulative energy preserved\n"
      "  -m pca     write/read pca information from this file\n"
      "  -n         normalize data before projection\n"
      "  -p idim    data input dimensions\n"
      "  -q odim    data output dimensions\n",
      prog, prog, prog, prog);
}

// input        -> (input) list of input file names
// block        -> (input) block size (number of rows to load in memory)
// exclude_dims -> (input) exclude these first/last dimensions from pca
// inp_dim      -> (input/output) number of input dimensions
// out_dim      -> (input/output) number of maximum output dimensions when
//                 doing projection
// miss_energy  -> (output) missed energy when projecting using all the
//                 selected eigenvectors
// eigval       -> (output) vector with the selected eigenvalues
//                 size: out_dim elements
// eigvec       -> (output) matrix with the selected eigenvectors from the data
//                 size: out_dim rows x (inp_dim - exclude_dims) columns
// mean         -> (output) vector with the mean of the input data
//                 size: inp_dim elements
// stddev       -> (output) vector with the standard deviation of the input data
//                 size: inp_dim elements
template <FORMAT_CODE fmt, typename real_t>
void compute_pca(
    vector<string> input, int block, int exclude_dims, int* inp_dim,
    int* out_dim, double* miss_energy, vector<real_t>* eigval,
    vector<real_t>* eigvec, vector<real_t>* mean, vector<real_t>* stddev) {
  int n = 0;  // number of data samples
  // process input to compute mean and co-moments
  compute_mean_comoments_from_inputs<fmt, real_t>(
      block, input, &n, inp_dim, mean, eigvec);
  CHECK_FMT(*inp_dim >= abs(exclude_dims),
            "Number of dimensions to exclude (%d) is bigger than the input "
            "dimensionality (%d)!", abs(exclude_dims), *inp_dim);
  CHECK_FMT(*inp_dim >= *out_dim,
            "Number of output dimensions (%d) is bigger than the input "
            "dimensionality (%d)!", *out_dim, *inp_dim);
  // if the number of output dimensions was not set, set it to the number of
  // input dimensions
  if (*out_dim <= 0) *out_dim = *inp_dim;
  // compute covariance from co-moments
  CHECK_FMT(n > 1, "You need at least 2 data points (only %d processed)!", n);
  for (int i = 0; i < (*inp_dim) * (*inp_dim); ++i) {
    (*eigvec)[i] /= (n - 1);
  }
  // compute standard deviation in each dimension
  stddev->resize(*inp_dim);
  for (int i = 0; i < (*inp_dim); ++i) {
    (*stddev)[i] = sqrt((*eigvec)[i * (*inp_dim) + i]);
  }
  // Up to here, we have the whole covariance matrix. But we will compute
  // only the eigenvalues and eigenvectors of the submatrix obtained after
  // removing the first/last `exclude_dims' rows and columns.
  const int pca_idim = *inp_dim - abs(exclude_dims);
  const int pca_odim = *out_dim - abs(exclude_dims);

  if (pca_odim > 0) {
    // Prepare memory for the eigenvalues
    eigval->resize(pca_idim);
    // Exclude from the covariance matrix the excluded dimensions.
    // NOTE: I am using the same space to store the covariance and the
    // eigenvalues, thus: this will destroy the covariance matrix!
    real_t* eigvec_offset = exclude_dims <= 0 ? eigvec->data() :
        eigvec->data() + exclude_dims * (*inp_dim) + exclude_dims;
    // Compute eigenvectors and eigenvalues
    CHECK(eig<real_t>(pca_idim, *inp_dim, eigvec_offset, eigval->data()) == 0);
    // Now, we have to computed the missed energy during projection, for the
    // dimensions that were not selected.
    *miss_energy = accumulate(eigval->begin() + pca_odim, eigval->end(), 0.0);
    eigval->resize(pca_odim);
    // Move all eigenvectors to the first rows, in order to free non-used
    // space of the covariance/eigenvectors matrix
    for (int r = 0; r < pca_odim; ++r) {
      for (int d = 0; d < pca_idim; ++d) {
        (*eigvec)[r * pca_idim + d] = eigvec_offset[r * (*inp_dim) + d];
      }
    }
    eigvec->resize(pca_odim * pca_idim);
  } else {
    eigval->resize(0);
    eigvec->resize(0);
    *miss_energy = pca_idim > 0 ? 1.0 : 0.0;
  }
}

template <FORMAT_CODE fmt, typename real_t>
void project_data(
    const vector<string>& input, const vector<string>& output,
    int block, int idim, int odim, int exclude_dims, double min_energy,
    bool normalize_data, const vector<real_t>& mean,
    const vector<real_t>& stddev, const vector<real_t>& eigval,
    const vector<real_t>& eigvec) {
  CHECK(input.size() > 0);
  CHECK(input.size() == output.size());
  CHECK(odim <= idim);
  CHECK(odim > 0 || min_energy > 0.0);
  CHECK_FMT(
      odim < 0 || odim >= exclude_dims,
      "Number of output dimensions (%d) is lower than the number of "
      "excluded dimensions from PCA (%d)!", odim, exclude_dims);
  const int max_odim_pca = eigval.size();
  CHECK(max_odim_pca > 0);
  // number of projected input and output dimensions
  int idim_pca = idim - abs(exclude_dims);
  int odim_pca = odim - abs(exclude_dims);
  // compute cumulative energy preserved and (optionally) output dimension
  const double total_energy = accumulate(eigval.begin(), eigval.end(), 0.0);
  // compute the number of output projected dimensions, if not computed yet,
  // and the preserved energy (cum_energy)
  double cum_energy = 0.0;
  if (odim_pca > 0) {
    // if the number of output dimensions is given, compute the preserved
    // energy
    cum_energy = accumulate(eigval.begin(), eigval.begin() + odim_pca, 0.0);
    cum_energy /= total_energy;
  } else if (odim <= 0) {
    // if the number of output dimensions is not given, then use as many
    // dimensions from pca as required to preserve the minimum required energy
    cum_energy = eigval[0];
    for (odim_pca = 1; odim_pca < max_odim_pca &&
             (cum_energy / total_energy < min_energy); ++odim_pca) {
      cum_energy += eigval[odim_pca];
    }
    odim = odim_pca + abs(exclude_dims);
    cum_energy = min(cum_energy / total_energy, 1.0);
  }
  // ----- process input files -----
  vector<real_t> x;  // data block
  vector<real_t> b;  // auxiliar data block
  if (idim > 0) {    // if we know the number of input dimensions, allocate mem
    x.resize(block * idim, 0);
    b.resize(block * idim, 0);
  }
  int n = 0;         // total number of processed samples (rows)
  for (size_t f = 0; f < input.size(); ++f) {
    // open input/output files
    const char* ifname = input[f] == "" ? "**stdin**" : input[f].c_str();
    const char* ofname = output[f] == "" ? "**stdout**" : output[f].c_str();
    FILE* ifile = input[f] == "" ? stdin : open_file(ifname, "rb");
    FILE* ofile = output[f] == "" ? stdout : open_file(ofname, "wb");
    // read input file header
    string h_name = "";
    int h_rows = -1, h_dims = -1;
    CHECK_FMT(
        !read_matrix_header<fmt>(ifile, &h_name, &h_rows, &h_dims),
        "Invalid header in file \"%s\"!", ifname);
    if (idim <= 0) {  // if we don't know yet the number of input dimensions
      CHECK_FMT(
          h_dims > 0, "Unknown number of dimensions in file \"%s\"!", ifname);
      idim = h_dims;  // the first input file determines the number of dims
      idim_pca = idim - abs(exclude_dims);
      CHECK(idim_pca <= odim_pca);
      x.resize(block * idim, 0);
      b.resize(block * idim, 0);
    } else {
      CHECK_FMT(
          idim == h_dims, "Bad number of dimensions in file \"%s\"!", ifname);
    }
    // write output header
    write_matrix_header<fmt>(ofile, h_name, h_rows, odim);
    // read, project and write data
    int fr = 0, be = 0, br = 0;
    while ((be = read_block<fmt, real_t>(ifile, block * idim, x.data())) > 0) {
      CHECK_FMT(be % idim == 0, "Corrupted matrix in file \"%s\"!", ifname);
      br = be / idim;
      fr += br;
      // project input data using pca
      project<real_t>(
          br, idim, odim, exclude_dims, eigvec.data(), mean.data(),
          normalize_data ? stddev.data() : NULL, x.data(), b.data());
      // output data
      write_matrix<fmt, real_t>(ofile, br, odim, x.data());
    }
    fclose(ifile);
    fclose(ofile);
    // update total number of processed rows
    n += fr;
    // if the number of read rows is not equal to the number of expected
    // rows, show a warning to the user
    if (h_rows > 0 && h_rows != fr) {
      WARN_FMT(
          "Number of processed rows (%d) is lower than expected (%d) "
          "in file \"%s\"!", fr, h_rows, ifname);
    }
  }
  // projection summary
  fprintf(stderr, "---------------- Projection summary -----------------\n");
  fprintf(stderr, "Processed rows: %d\n", n);
  fprintf(stderr, "Input dimension: %d\n", idim);
  fprintf(stderr, "Output dimension: %d\n", odim);
  if (exclude_dims <= 0) {
    fprintf(
        stderr, "Projected dimensions: %d-%d\n", abs(exclude_dims) + 1, idim);
  } else {
    fprintf(
        stderr, "Projected dimensions: %d-%d\n", 1, idim_pca);
  }
  fprintf(stderr, "Preserved energy: %.4g%%\n", cum_energy * 100.0);
  fprintf(stderr, "-----------------------------------------------------\n");
}

template <typename real_t>
void pca_summary(
    int inp_dim, int out_dim, int exclude_dims, double miss_energy,
    const vector<real_t>& eigval) {
  const int pca_idims = inp_dim - abs(exclude_dims);
  const double total_energy = accumulate(
      eigval.begin(), eigval.end(), miss_energy);
  // compute energy quantiles
  // (how many dimensions you need to preserve % of energy)
  const double quant_val[] = {0.25, 0.5, 0.75, 1.0};
  int quant_dim[] = {pca_idims, pca_idims, pca_idims, pca_idims};
  double cum_energy = 0.0;
  for (int i = 0, q = 0; i < pca_idims && q < 4; ++i) {
    cum_energy += eigval[i];
    for (; q < 4 && cum_energy >= total_energy * quant_val[q]; ++q) {
      quant_dim[q] = i + 1;
    }
  }
  fprintf(
      stderr,
      "-------------------- PCA summary --------------------\n"
      "Input dimensions: %d\n"
      "Non-projected dimensions: %d\n"
      "Energy quantiles: 25%% -> %d, 50%% -> %d, 75%% -> %d, 100%% -> %d\n"
      "-----------------------------------------------------\n",
      inp_dim, exclude_dims, quant_dim[0], quant_dim[1], quant_dim[2],
      quant_dim[3]);
}

template <FORMAT_CODE fmt, typename real_t>
void do_work(
    const bool do_compute_pca, const bool do_project_data, const string& pca_fn,
    const vector<string>& input, const vector<string>& output, int block,
    int inp_dim, int out_dim, double min_energy, bool normalize_data,
    int exclude_dims) {
  vector<real_t> mean;
  vector<real_t> stdev;
  vector<real_t> eigval;
  vector<real_t> eigvec;
  double miss_energy = 0.0;
  if (do_compute_pca) {
    compute_pca<fmt, real_t>(
        input, block, exclude_dims, &inp_dim, &out_dim, &miss_energy,
        &eigval, &eigvec, &mean, &stdev);
    if (!do_project_data || pca_fn != "") {
      save_pca<real_t>(
          pca_fn, exclude_dims, miss_energy, mean, stdev, eigval, eigvec);
    }
  } else {
    CHECK_MSG(pca_fn != "", "Specify a pca file to load from!");
    if (exclude_dims != 0) {
      WARN_FMT(
          "Ignoring \"-e %d\": excluded dimensions will be read from the "
          "pca file...", exclude_dims);
    }
    load_pca<real_t>(
        pca_fn, &inp_dim, &out_dim, &exclude_dims, &miss_energy,
        &mean, &stdev, &eigval, &eigvec);
  }
  pca_summary(inp_dim, out_dim, exclude_dims, miss_energy, eigval);
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

  double min_energy = -1.0;
  bool do_compute_pca = false;
  bool do_project_data = false;
  FORMAT_CODE format = FMT_ASCII;
  const char* format_str = NULL;
  while ((opt = getopt(argc, argv, "CPb:de:f:hj:m:np:q:")) != -1) {
    switch (opt) {
      case 'C':
        do_compute_pca = true;
        break;
      case 'P':
        do_project_data = true;
        break;
      case 'b':
        block = atoi(optarg);
        CHECK_FMT(block > 0, "Block size must be positive (-b %d)!", block);
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
        CHECK_FMT(format != FMT_UNKNOWN, "Unknown format (-f \"%s\")!", optarg);
        break;
      case 'h':
        help(argv[0]);
        return 0;
      case 'j':
        min_energy = atof(optarg);
        CHECK_FMT(
            min_energy >= 0.0 && min_energy <= 1.0,
            "Invalid minimum cumulative energy (-j %f)!", min_energy);
        break;
      case 'm':
        pca_fn = optarg;
        break;
      case 'p':
        inp_dim = atoi(optarg);
        CHECK_FMT(
            inp_dim > 0, "Input dimension must be positive (-p %d)!", inp_dim);
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
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  if (min_energy < 0) min_energy = 1.0;

  // input & output file names
  vector<string> input, output;
  if (do_project_data) {
    for (int a = optind; a < argc; a+=2) {
      input.push_back(argv[a]);
      if (a + 1 < argc) {
        output.push_back(argv[a + 1]);
      } else {
        output.push_back("");
      }
    }
  } else {
    for (int a = optind; a < argc; ++a) {
      input.push_back(argv[a]);
    }
  }
  // when reading from stdin, we cannot process data twice!
  CHECK_MSG(
      !do_compute_pca || !do_project_data || input.size() > 0,
      "You cannot perform PCA and project data when reading from stdin!");

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
        do_work<FMT_BINARY, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      } else {
        do_work<FMT_BINARY, double>(
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
    case FMT_PRHLT_HTK:
      if (simple_precision) {
        do_work<FMT_PRHLT_HTK, float>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_energy, normalize_data,
            exclude_dims);
      } else {
        do_work<FMT_PRHLT_HTK, double>(
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
