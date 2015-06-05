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
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/pca.h"
#include "fast_pca/fast_pca_common.h"
#include "fast_pca/logging.h"

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
      "  -e dims    do not project first (positive) or last (negative) dims\n"
      "  -f format  format of the data matrix (ascii, binary, octave, vbosch,\n"
      "             htk, mat4)\n"
      "  -j energy  minimum relative amount of energy preserved\n"
      "  -m pca     write/read pca information to/from this file\n"
      "  -n         normalize data before projection\n"
      "  -p idim    data input dimensions\n"
      "  -q odim    data output dimensions\n",
      prog, prog, prog, prog);
}

// input          -> (input) list of input file names
// block          -> (input) block size (number of rows to load in memory)
// exclude_dims   -> (input) exclude these first/last dimensions from pca
// min_rel_energy -> (input) minimum amount of relative energy to preserve
// inp_dim        -> (input/output) number of input dimensions
// out_dim        -> (input/output) number of maximum output dimensions when
//                   doing projection
// miss_energy    -> (output) missed energy when projecting using all the
//                   selected eigenvectors
// eigval         -> (output) vector with the selected eigenvalues
//                   size: out_dim elements
// eigvec         -> (output) matrix with the selected eigenvectors
//                   size: out_dim rows x (inp_dim - exclude_dims) columns
// mean           -> (output) vector with the mean
//                   size: inp_dim elements
// stddev         -> (output) vector with the standard deviation
//                   size: inp_dim elements
template <FORMAT_CODE fmt, typename real_t>
void compute_pca(
    const vector<string>& input, int block, int exclude_dims,
    double min_rel_energy, int* inp_dim, int* out_dim, double* miss_energy,
    vector<real_t>* eigval, vector<real_t>* eigvec, vector<real_t>* mean,
    vector<real_t>* stddev) {
  int n = 0;  // number of data samples
  // process input to compute mean and co-moments
  compute_mean_comoments_from_inputs<fmt, real_t>(
      block, input, &n, inp_dim, mean, eigvec);
  CHECK_FMT(*inp_dim >= *out_dim,
            "Number of output dimensions (%d) is bigger than the input "
            "dimensions (%d)!", *out_dim, *inp_dim);
  CHECK_FMT(*inp_dim >= abs(exclude_dims),
            "Number of non-projected dimensions (%d) is bigger than the input "
            "dimensions (%d)!", abs(exclude_dims), *inp_dim);
  CHECK_FMT(*out_dim < 1 || *out_dim >= abs(exclude_dims),
            "Number of non-projected dimensions (%d) is bigger than the output "
            "dimensions (%d)!", abs(exclude_dims), *out_dim);
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
  // compute eigenvectors and eigenvalues of the covariance matrix
  compute_pca_from_covariance<real_t>(
      exclude_dims, min_rel_energy, *inp_dim, out_dim, miss_energy,
      eigvec, eigval);
}

template <FORMAT_CODE fmt, typename real_t>
int project_data(
    const vector<string>& input, const vector<string>& output,
    const int block, const int odim, const int exclude_dims,
    const bool normalize_data, const vector<real_t>& mean,
    const vector<real_t>& stddev, const vector<real_t>& eigval,
    const vector<real_t>& eigvec) {
  const int idim = mean.size();
  CHECK(idim > 0);
  CHECK(odim > 0);
  CHECK(odim <= idim);
  CHECK(input.size() > 0);
  CHECK(input.size() == output.size());
  // ----- process input files -----
  vector<real_t> x(block * idim, 0);  // data block
  vector<real_t> z(block * odim, 0);  // auxiliar data block
  // matrix reader
  unique_ptr<MatrixFile> mr(MatrixFile::Create<fmt>());  // matrix reader
  unique_ptr<MatrixFile> mw(MatrixFile::Create<fmt>());  // matrix writer

  int n = 0;         // total number of processed samples (rows)
  for (size_t f = 0; f < input.size(); ++f) {
    // open input/output files
    const char* ifname = input[f] == "" ? "**stdin**" : input[f].c_str();
    const char* ofname = output[f] == "" ? "**stdout**" : output[f].c_str();
    FILE* ifile = input[f] == "" ? stdin : open_file(ifname, "rb");
    FILE* ofile = output[f] == "" ? stdout : open_file(ofname, "wb");
    // read input file header
    mr->file(ifile);
    CHECK_FMT(mr->read_header(), "Invalid header in file \"%s\"!", ifname);
    CHECK_FMT(
        mr->cols() < 0 || mr->cols() == idim,
        "Bad number of dimensions in file \"%s\" (found: %d, expected: %d)!",
        ifname, mr->cols(), idim);
    // write output file header
    mw->file(ofile);
    mw->copy_header_from(*mr);
    mw->cols(odim);
    mw->write_header();
    // read, project and write data
    int fr = 0, be = 0, br = 0;
    while ((be = mr->read_block(block * idim, x.data())) > 0) {
      CHECK_FMT(
          be % idim == 0,
          "Corrupted matrix in file \"%s\" (block expected a multiple of "
          "%d elements, but %d where read)!\n", ifname, idim, be);
      br = be / idim;
      fr += br;
      // project input data using pca
      project<real_t>(
          br, idim, odim, exclude_dims, eigvec.data(), mean.data(),
          normalize_data ? stddev.data() : NULL, x.data(), z.data());
      // output data
      mw->write_block(br * odim, z.data());
    }
    fclose(ifile);
    fclose(ofile);
    // update total number of processed rows
    n += fr;
    // if the number of read rows is not equal to the number of expected
    // rows, show a warning to the user
    if (mr->rows() > 0 && mr->rows() != fr) {
      WARN_FMT(
          "Number of processed rows (%d) is lower than expected (%d) "
          "in file \"%s\"!", fr, mr->rows(), ifname);
    }
  }
  return n;
}

void projection_summary(
    const int n, const int idim, const int odim, const int exclude_dims,
    const double miss_energy, const double kept_energy) {
  const double rel_kept_energy = kept_energy / (miss_energy + kept_energy);
  fprintf(stderr, "---------------- Projection summary -----------------\n");
  fprintf(stderr, "Processed rows: %d\n", n);
  fprintf(stderr, "Input dimension: %d\n", idim);
  fprintf(stderr, "Output dimension: %d\n", odim);
  if (exclude_dims > 0) {
    fprintf(
        stderr, "Projected dimensions: %d-%d\n", 1 + exclude_dims, idim);
  } else {
    fprintf(
        stderr, "Projected dimensions: %d-%d\n", 1, idim + exclude_dims);
  }
  fprintf(stderr, "Preserved energy: %.4g%%\n", rel_kept_energy * 100.0);
  fprintf(stderr, "-----------------------------------------------------\n");
}

template <FORMAT_CODE fmt, typename real_t>
void do_work(
    const bool do_compute_pca, const bool do_project_data,
    const string& pca_fn,
    const vector<string>& input, const vector<string>& output, int block,
    int inp_dim, int out_dim, double min_rel_energy, bool normalize_data,
    int exclude_dims) {
  vector<real_t> mean;
  vector<real_t> stdev;
  vector<real_t> eigval;
  vector<real_t> eigvec;
  double miss_energy = 0.0;
  if (do_compute_pca) {
    // Compute PCA from input files
    compute_pca<fmt, real_t>(
        input, block, exclude_dims, min_rel_energy, &inp_dim,
        &out_dim, &miss_energy, &eigval, &eigvec, &mean, &stdev);
    if (!do_project_data || pca_fn != "") {
      save_pca<real_t>(
          pca_fn, exclude_dims, miss_energy, mean, stdev, eigval, eigvec);
    }
  } else {
    // Load PCA projection info from a previously generated file
    CHECK_MSG(pca_fn != "", "Specify a pca file to load from!");
    if (exclude_dims != 0) {
      WARN_FMT(
          "Ignoring \"-e %d\": non-projected dimensions will be read from the "
          "pca file...", exclude_dims);
    }
    load_pca<real_t>(
        pca_fn, &exclude_dims, &miss_energy, &mean, &stdev, &eigval, &eigvec);
    CHECK_FMT(
        inp_dim < 1 || inp_dim == mean.size(),
        "Number of input dimensions (%d) does not match to the number of "
        "dimensions read from the pca file (%d)!", inp_dim, (int)mean.size());
    inp_dim = inp_dim > 0 ? inp_dim : mean.size();
  }
  // Compute PCA summary
  vector<real_t> cumulative_energy;
  compute_cumulative_energy(eigval, &cumulative_energy);
  pca_summary<real_t>(
      inp_dim, exclude_dims, miss_energy, cumulative_energy);
  if (do_project_data) {
    // Project input data using the PCA information
    const double total_energy = miss_energy + cumulative_energy.back();
    int pca_odim = 0;
    if (out_dim < 1 && min_rel_energy > 0.0) {
      pca_odim = compute_pca_output_dim<real_t>(
          cumulative_energy, min_rel_energy, miss_energy);
      out_dim = pca_odim + abs(exclude_dims);
    } else if (out_dim > 0) {
      pca_odim = out_dim - abs(exclude_dims);
    } else {
      out_dim = inp_dim;
      pca_odim = out_dim - abs(exclude_dims);
    }
    miss_energy = total_energy - cumulative_energy[pca_odim];
    const int n = project_data<fmt, real_t>(
        input, output, block, out_dim, exclude_dims, normalize_data, mean,
        stdev, eigval, eigvec);
    projection_summary(
        n, inp_dim, out_dim, exclude_dims, miss_energy,
        cumulative_energy[pca_odim]);
  }
}

int main(int argc, char** argv) {
  int opt = -1;
  int inp_dim = -1, out_dim = -1;
  int exclude_dims = 0;
  int block = 1000;
  bool simple_precision = true;
  bool normalize_data = false;
  bool do_compute_pca = false;
  bool do_project_data = false;
  double min_rel_energy = -1.0;
  string pca_fn = "";
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
        min_rel_energy = atof(optarg);
        CHECK_FMT(
            min_rel_energy >= 0.0 && min_rel_energy <= 1.0,
            "Invalid minimum amount of relative energy (-j %f)!",
            min_rel_energy);
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
  if (format_str) fprintf(stderr, " -f \"%s\"", format_str);
  if (min_rel_energy > 0) fprintf(stderr, " -j %g", min_rel_energy);
  if (pca_fn != "") fprintf(stderr, " -m \"%s\"", pca_fn.c_str());
  if (normalize_data) fprintf(stderr, " -n");
  if (inp_dim > 0) fprintf(stderr, " -p %d", inp_dim);
  if (out_dim > 0) fprintf(stderr, " -q %d", out_dim);
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

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
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      } else {
        do_work<FMT_ASCII, double>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      }
      break;
    case FMT_BINARY:
      if (simple_precision) {
        do_work<FMT_BINARY, float>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      } else {
        do_work<FMT_BINARY, double>(
            do_compute_pca, do_project_data, pca_fn, input,
            output, block, inp_dim, out_dim, min_rel_energy, normalize_data,
            exclude_dims);
      }
      break;
    case FMT_OCTAVE:
      if (simple_precision) {
        do_work<FMT_OCTAVE, float>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      } else {
        do_work<FMT_OCTAVE, double>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      }
      break;
    case FMT_VBOSCH:
      if (simple_precision) {
        do_work<FMT_VBOSCH, float>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      } else {
        do_work<FMT_VBOSCH, double>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      }
      break;
    case FMT_HTK:
      if (simple_precision) {
        do_work<FMT_HTK, float>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      } else {
        do_work<FMT_HTK, double>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      }
      break;
    case FMT_MAT4:
      if (simple_precision) {
        do_work<FMT_MAT4, float>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      } else {
        do_work<FMT_MAT4, double>(
            do_compute_pca, do_project_data, pca_fn, input, output, block,
            inp_dim, out_dim, min_rel_energy, normalize_data, exclude_dims);
      }
      break;
    default:
      ERROR("Not implemented for this format!");
  }
  return 0;
}
