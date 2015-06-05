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

#ifndef FAST_PCA_FAST_PCA_COMMON_H_
#define FAST_PCA_FAST_PCA_COMMON_H_

#include "fast_pca/file.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/math.h"
#include "fast_pca/pca.h"

#include <memory>
#include <string>
#include <vector>

using std::string;
using std::unique_ptr;
using std::vector;

template <FORMAT_CODE fmt, typename real_t>
void compute_mean_comoments_from_inputs(
    int block, vector<string> input, int* n, int* inp_dim,
    vector<real_t>* M, vector<real_t>* C) {
  CHECK(!input.empty());
  CHECK(block > 0);
  const vector<real_t> ones(block, 1);
  vector<real_t> x;  // data block
  vector<real_t> m;  // mean of the current block
  vector<real_t> d;  // diff between global and block mean
  if (*inp_dim > 0) {
    x.resize(block * *inp_dim, 0);
    m.resize(*inp_dim, 0);
    d.resize(*inp_dim, 0);
    M->resize(*inp_dim, 0);
    C->resize((*inp_dim) * (*inp_dim), 0);
  }
  *n = 0;            // total processed rows
  unique_ptr<MatrixFile> mh(MatrixFile::Create<fmt>());
  for (size_t f = 0; f < input.size(); ++f) {
    const char* fname = input[f] == "" ? "**stdin**" : input[f].c_str();
    FILE* file = input[f] == "" ? stdin : open_file(fname, "rb");
    mh->file(file);
    CHECK_FMT(
        mh->read_header(), "Failed to read header in file \"%s\"!", fname);
    if (*inp_dim < 1) {
      CHECK_FMT(
          mh->cols() > 0,
          "Number of input dimensions could not be determined by file \"%s\" "
          "(number of read columns in file: %d)!", fname, mh->cols());
      *inp_dim = mh->cols();
      x.resize(block * (*inp_dim), 0);
      m.resize(*inp_dim, 0);
      d.resize(*inp_dim, 0);
      M->resize(*inp_dim, 0);
      C->resize((*inp_dim) * (*inp_dim), 0);
    } else {
      CHECK_FMT(
          mh->cols() < 0 || mh->cols() == *inp_dim,
          "Number of read dimensions in file \"%s\" (%d) is not the "
          "expected (%d)!", fname, mh->cols(), *inp_dim);
    }
    int fr = 0, be = 0, br = 0;
    while ((be = mh->read_block(block * (*inp_dim), x.data())) > 0) {
      CHECK_FMT(
        be % (*inp_dim) == 0,
        "Corrupted matrix in file \"%s\" (block expected a multiple of "
        "%d elements, but %d where read)!\n",
        fname, *inp_dim, be);
      br = be / (*inp_dim);
      fr += br;
      // compute block mean
      gemv<real_t>(
          'T', br, *inp_dim, 1.0 / br, x.data(), *inp_dim, ones.data(), 1,
          0, m.data(), 1);
      // subtract mean to the current block
      for (int i = 0; i < br; ++i) {
        axpy<real_t>(*inp_dim, -1, m.data(), x.data() + i * (*inp_dim));
      }
      // d = M - m
      memcpy(d.data(), M->data(), sizeof(real_t) * (*inp_dim));
      axpy<real_t>(*inp_dim, -1, m.data(), d.data());
      // update co-moments matrix
      // C += (x - m)' * (x - m)
      gemm<real_t>(
          'T', 'N', *inp_dim, *inp_dim, br, 1, x.data(), *inp_dim, x.data(),
          *inp_dim, 1, C->data(), *inp_dim);
      // C += D * D' * (br * n) / (br + n)
      const int nn = *n + br;
      const real_t cf = br * ((*n) / (1.0 * nn));
      ger<real_t>(*inp_dim, *inp_dim, cf, d.data(), d.data(), C->data());
      // update mean
      for (int i = 0; i < *inp_dim; ++i) {
        (*M)[i] = ((*n) * (*M)[i] + br * m[i]) / nn;
      }
      // update total number of processed rows
      *n = nn;
    }
    fclose(file);
  }
}


template <typename real_t>
void compute_cumulative_energy(
    const vector<real_t>& eigval, vector<real_t>* cumulative_energy) {
  cumulative_energy->resize(eigval.size() + 1);
  (*cumulative_energy)[0] = 0.0;
  for (size_t k = 1; k <= eigval.size(); ++k) {
    (*cumulative_energy)[k] = (*cumulative_energy)[k - 1] + \
        (eigval[k] > 0.0 ? fabs(eigval[k]) : 0.0);
  }
}


template <typename real_t>
int compute_pca_output_dim(
    const vector<real_t>& cumulative_energy, const double min_rel_energy,
    const double miss_energy) {
  CHECK(min_rel_energy > 0.0);
  size_t k = 0;
  const double total_energy = miss_energy + cumulative_energy.back();
  for (; k + 1 < cumulative_energy.size() &&
          cumulative_energy[k] < min_rel_energy * total_energy; ++k) {}
  return k;
}


template <typename real_t>
void count_negative_and_zero_eigenvalues(
    const vector<real_t>& eigval, int* num_neg, int* num_zero) {
  *num_neg = 0;
  *num_zero = 0;
  size_t k = eigval.size();
  for (; k > 0 && eigval[k - 1] < 0.0; --k, ++(*num_neg)) {}
  for (; k > 0 && eigval[k - 1] == 0.0; --k, ++(*num_zero)) {}
}


template <typename real_t>
void compute_pca_from_covariance(
    const int exclude_dims, const double min_rel_energy, const int inp_dim,
    int* out_dim, double* miss_energy, vector<real_t>* eigvec,
    vector<real_t>* eigval) {
  const int pca_idim = inp_dim - abs(exclude_dims);
  if (pca_idim < 1) {
    WARN_FMT(
        "All input dimensions (%d) were excluded from pca. "
        "The output will just be the mean centered input data.", inp_dim);
    eigvec->resize(0);
    eigval->resize(0);
    *miss_energy = 0.0;
    *out_dim = *out_dim > 0 ? *out_dim : inp_dim;
  } else {
    // Exclude from the non-projected dimensions.
    real_t* eigvec_ptr = eigvec->data() + \
        (exclude_dims > 0) * exclude_dims * (inp_dim + 1);
    // Prepare memory for the eigenvalues
    eigval->resize(pca_idim);
    // Compute eigenvectors and eigenvalues
    CHECK(eig<real_t>(pca_idim, inp_dim, eigvec_ptr, eigval->data()) == 0);
    // Check for zero or negative eigenvalues
    int num_neg_eigval = 0, num_zero_eigval = 0;
    count_negative_and_zero_eigenvalues<real_t>(
        *eigval, &num_neg_eigval, &num_zero_eigval);
    if (num_zero_eigval > 0 || num_neg_eigval > 0) {
      WARN_FMT(
          "Covariance matrix is not positive definite: %d zero and %d negative "
          "eigenvalues found! This probably means there is an strong "
          "correlation between some of your input variables.",
          num_zero_eigval, num_neg_eigval);
      if (-eigval->back() > 1E-6) {
        WARN_FMT(
            "The lowest eigenvalue (%g) is far from zero. This may be a "
            "serius problem caused by numerical precision errors. Try using "
            "the `-d' option for higher precision computations.",
            eigval->back());
      }
    }
    // compute cumulative energy achieved by adding each eigenvector
    vector<real_t> cumulative_energy;
    compute_cumulative_energy(*eigval, &cumulative_energy);
    const double total_energy = cumulative_energy.back();
    // compute pca output dimensions, if not given
    int pca_odim = *out_dim - abs(exclude_dims);
    if (*out_dim < 1 && min_rel_energy > 0.0) {
      // number of output dimensions was not given, but minimum
      // relative energy to preserve was given instead
      pca_odim = compute_pca_output_dim<real_t>(
          cumulative_energy, min_rel_energy, 0.0);
      *out_dim = pca_odim + abs(exclude_dims);
    } else if (*out_dim < 1) {
      // neither number of output dimensions or minimum relative energy
      // to preserve was given
      pca_odim = pca_idim;
      *out_dim = pca_odim + abs(exclude_dims);
    }
    // missed energy during pca projection
    *miss_energy = total_energy - cumulative_energy[pca_odim];
    eigval->resize(pca_odim);
    // Move all eigenvectors to the first rows, in order to free non-used
    // space of the eigenvectors matrix
    for (int r = 0; r < pca_odim; ++r) {
      for (int d = 0; d < pca_idim; ++d) {
        (*eigvec)[r * pca_idim + d] = eigvec_ptr[r * inp_dim + d];
      }
    }
    eigvec->resize(pca_odim * pca_idim);
  }
}


template <typename real_t>
void pca_summary(
    const int inp_dim, const int exclude_dims, const double miss_energy,
    const vector<real_t>& cumulative_energy) {
  if (cumulative_energy.size() > 1) {
    const double total_energy = cumulative_energy.back() + miss_energy;
    const double rel_kept_energy = cumulative_energy.back() / total_energy;
    const int max_pca_odim = cumulative_energy.size() - 1;
    const double quant_val[] = {0.25, 0.5, 0.75, 1.0};
    int quant_dim[] = {max_pca_odim, max_pca_odim, max_pca_odim, max_pca_odim};
    for (int i = 1, q = 0; i <= max_pca_odim && q < 4; ++i) {
      if (cumulative_energy[i] < quant_val[q] * total_energy) continue;
      quant_dim[q++] = i;
    }
    fprintf(
        stderr,
        "-------------------- PCA summary --------------------\n"
        "Input dimensions: %d\n"
        "Projectable input dimensions: %d-%d\n"
        "Maximum output dimensions: %d\n"
        "Maximum preservervable energy: %.4g%%\n"
        "Energy quantiles: %d%% -> %d, %d%% -> %d, %d%% -> %d, %d%% -> %d\n"
        "-----------------------------------------------------\n",
        inp_dim,
        exclude_dims < 0 ? 1 : exclude_dims + 1,
        exclude_dims < 0 ? inp_dim + exclude_dims : inp_dim,
        max_pca_odim, rel_kept_energy * 100.0,
        static_cast<int>(round(quant_val[0] * rel_kept_energy * 100.0)),
        quant_dim[0],
        static_cast<int>(round(quant_val[1] * rel_kept_energy * 100.0)),
        quant_dim[1],
        static_cast<int>(round(quant_val[2] * rel_kept_energy * 100.0)),
        quant_dim[2],
        static_cast<int>(round(quant_val[3] * rel_kept_energy * 100.0)),
        quant_dim[3]);
  } else {
    fprintf(
        stderr,
        "-------------------- PCA summary --------------------\n"
        "Input dimensions: %d\n"
        "Projectable input dimensions: None\n"
        "-----------------------------------------------------\n",
        inp_dim);
  }
}

#endif  // FAST_PCA_FAST_PCA_COMMON_H_
