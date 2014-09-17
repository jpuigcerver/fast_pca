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
#include <cstring>
#include <string>
#include <vector>

#include "fast_pca/file.h"
#include "fast_pca/file_pca.h"
#include "fast_pca/math.h"

using std::string;
using std::vector;

void help(const char* prog) {
  fprintf(
      stderr,
      "Usage: %s [-b size] [-d] [-f format] [-o output] [-d dim] [input ...]\n"
      "Options:\n"
      "  -b size    process this number of rows at once\n"
      "  -d         use double precision\n"
      "  -f format  format of the data matrix\n"
      "  -o output  output file\n"
      "  -p dim     data dimensions\n",
      prog);
}

template <FORMAT_CODE fmt, typename real_t>
void do_work(int block, int dims, string output, vector<string>* input) {
  // open input files
  vector<FILE*> files;
  open_files("rb", "**stdin**", stdin, input, &files);
  // process files
  const vector<real_t> ones(block, 1);
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = (*input)[f].c_str();
    FILE* file = files[f];
    int h_rows = -1, h_dims = -1;
    if (read_matrix_header<fmt>(file, NULL, &h_rows, &h_dims)) {
      fprintf(stderr, "ERROR: Invalid header in file \"%s\"!\n", fname);
      exit(1);
    }
    if (dims < 0) dims = h_dims;
    else if (dims != h_dims) {
      fprintf(
          stderr, "ERROR: Bad number of dimensions in file \"%s\"!\n", fname);
      exit(1);
    }
  }
  vector<real_t> x(block * dims, 0);  // data block
  vector<real_t> m(dims, 0);          // mean of the current block
  vector<real_t> D(dims, 0);          // diff between global and block mean
  vector<real_t> M(dims, 0);          // global mean
  vector<real_t> C(dims * dims, 0);   // global co-moments matrix
  int n = 0;                          // total processed rows
  for (size_t f = 0; f < files.size(); ++f) {
    const char* fname = (*input)[f].c_str();
    FILE* file = files[f];
    int fr = 0, be = 0, br = 0;
    while((be = read_block<fmt, real_t>(file, block * dims, x.data())) > 0) {
      if (be % dims != 0) {
        fprintf(stderr, "ERROR: Corrupted matrix in file \"%s\"!\n", fname);
        exit(1);
      }
      br = be / dims;
      fr += br;

      // compute block mean
      gemv<real_t>(
          'T', br, dims, 1.0 / br, x.data(), dims, ones.data(), 1,
          0, m.data(), 1);
      // subtract mean to the current block
      for (int i = 0; i < br; ++i) {
        axpy<real_t>(dims, -1, m.data(), x.data() + i * dims);
      }
      // D = M - m
      memcpy(D.data(), M.data(), sizeof(real_t) * dims);
      axpy<real_t>(dims, -1, m.data(), D.data());
      // update co-moments matrix
      // C += (x - m)' * (x - m)
      gemm<real_t>(
          'T', 'N', dims, dims, br, 1, x.data(), dims, x.data(), dims, 1.0,
          C.data(), dims);
      // C += D * D' * (br * n) / (br + n)
      ger<real_t>(
          dims, dims, (1.0 * br) * n / (n + br), D.data(), D.data(), C.data());
      // update mean
      for (int d = 0; d < dims; ++d) {
        M[d] = (n * M[d] + br * m[d]) / (n + br);
      }
      // update total number of processed rows
      n += br;
    }
  }
  // close input files
  close_files(files);
  // output number of processed rows, mean and co-moments matrix
  save_n_mean_cov(output, n, dims, M, C);
}

int main(int argc, char** argv) {
  int opt = -1;
  int dims = -1;             // number of dimensions
  int block = 1000;          // block size
  bool simple = true;        // use simple precision ?
  string output = "";
  FORMAT_CODE format = FMT_VBOSCH;
  const char* format_str = NULL;

  while ((opt = getopt(argc, argv, "db:f:o:p:h")) != -1) {
    switch (opt) {
      case 'd':
        simple = false;
        break;
      case 'b':
        block = atoi(optarg);
        if (block < 1) {
          fprintf(stderr, "ERROR: Block size must be positive!\n");
          exit(1);
        }
        break;
      case 'f':
        format_str = optarg;
        format = format_code_from_name(format_str);
        if (format == FMT_UNKNOWN) {
          fprintf(stderr, "ERROR: Unknown format \"%s\"!\n", optarg);
          exit(1);
        }
        break;
      case 'o':
        output = optarg;
        break;
      case 'p':
        dims = atoi(optarg);
        if (dims < 1) {
          fprintf(stderr, "ERROR: Input dimension must be positive!\n");
          exit(1);
        }
        break;
      case 'h':
        help(argv[0]);
        return 0;
      default:
        return 1;
    }
  }

  fprintf(stderr, "-------------------- Command line -------------------\n");
  fprintf(stderr, "%s", argv[0]);
  fprintf(stderr, " -b %d", block);
  if (!simple) fprintf(stderr, " -d");
  if (format_str) fprintf(stderr, " -f \"%s\"", format_str);
  if (output != "") fprintf(stderr, "-o %s", output.c_str());
  if (dims > 0) fprintf(stderr, " -p %d", dims);
  for (int a = optind; a < argc; ++a) {
    fprintf(stderr, " \"%s\"", argv[a]);
  }
  fprintf(stderr, "\n-----------------------------------------------------\n");

  vector<string> input;
  for (int a = optind; a < argc; ++a) { input.push_back(argv[a]); }

  switch (format) {
    case FMT_ASCII:
      if (simple) do_work<FMT_ASCII, float>(block, dims, output, &input);
      else do_work<FMT_ASCII, double>(block, dims, output, &input);
      break;
    case FMT_BINARY:
      if (simple) do_work<FMT_BINARY, float>(block, dims, output, &input);
      else do_work<FMT_BINARY, double>(block, dims, output, &input);
      break;
    case FMT_OCTAVE:
      if (simple) do_work<FMT_OCTAVE, float>(block, dims, output, &input);
      else do_work<FMT_OCTAVE, double>(block, dims, output, &input);
      break;
    case FMT_VBOSCH:
      if (simple) do_work<FMT_VBOSCH, float>(block, dims, output, &input);
      else do_work<FMT_VBOSCH, double>(block, dims, output, &input);
      break;
    default:
      fprintf(stderr, "ERROR: Not implemented for this format!\n");
      exit(1);
  }

  return 0;
}
