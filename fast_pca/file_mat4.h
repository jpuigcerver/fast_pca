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

#ifndef FAST_PCA_FILE_MAT4_H_
#define FAST_PCA_FILE_MAT4_H_

#include "fast_pca/file.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

class MatrixFile_MAT4 : public MatrixFile {
 protected:
  string name_;
  uint32_t mopt_;
  uint32_t prec_;
  char order_;
  bool swap_;

  // this is used to map a given type T to its precision ID,
  // according to MAT-v4 format file.
  // Example: float_precision = type_prec<float>::prec;
  template <typename T>
  struct type2prec { static uint8_t prec; };

 public:
  MatrixFile_MAT4() : MatrixFile(FMT_MAT4) {}
  explicit MatrixFile_MAT4(FILE* file) :
      MatrixFile(file), mopt_(0), prec_(0), order_(0), swap_(0) {}
  explicit MatrixFile_MAT4(
      FILE* file, int rows, int cols, const string& name, uint8_t prec) :
      MatrixFile(file), name_(name), mopt_(0), prec_(prec), order_(0),
      swap_(0) {
    MatrixFile::rows_ = rows;
    MatrixFile::cols_ = cols;
    CHECK(prec_ < 10);
  }

  inline void name(const string& name) { name_ = name; }
  inline const string& name() const { return name_; }

  virtual bool copy_header_from(const MatrixFile& other);
  virtual bool read_header();
  virtual void write_header() const;

  // Trick to avoid duplication of code all read/write methods for
  // different types are very similar. However virtual templated functions
  // are not allowed. To avoid this, I create a non-virtual templated function
  // with the same name, then the virtual functions call their respective
  // templated functions.
  template <typename T>
  int read_block(int n, T* m) const;
  template <typename T>
  void write_block(int n, const T* m) const;

  virtual int read_block(int n, float* m) const {
    return read_block<float>(n, m);
  }
  virtual int read_block(int n, double* m) const {
    return read_block<double>(n, m);
  }
  virtual void write_block(int n, const float* m) const {
    write_block<float>(n, m);
  }
  virtual void write_block(int n, const double* m) const {
    write_block<double>(n, m);
  }

  // load/save complete scalar/matrix data
  template <typename T>
  static void load(FILE* file, string* name, T* n) {
    MatrixFile_MAT4 h(file);
    // read and check header
    CHECK(h.read_header());
    CHECK_FMT(
        h.rows_ == 1 && h.cols_ == 1,
        "Matrix \"%s\" size (%dx%d) is not the expected one (1x1)",
        h.name_.c_str(), h.rows_, h.cols_);
    // read data
    CHECK(h.read_block(1, n) == 1);
    *name = h.name_;
  }

  template <typename T>
  static void load(
      FILE* file, string* name, int* rows, int* cols, vector<T>* m) {
    MatrixFile_MAT4 h(file);
    h.read_header();
    CHECK(h.rows_ >= 0 && h.cols_ >= 0);
    m->resize(h.rows_ * h.cols_);
    CHECK(h.read_block(h.rows_ * h.cols_, m->data()) == h.rows_ * h.cols_);
    *rows = h.rows_;
    *cols = h.cols_;
    *name = h.name_;
  }

  template <typename T>
  static void save(FILE* file, const string& name, T n) {
    MatrixFile_MAT4 h(file, 1, 1, name, type2prec<T>::prec);
    h.write_header();
    h.write_block(1, &n);
  }

  template <typename T>
  static void save(
      FILE* file, const string& name, int rows, int cols,
      const vector<T>& m) {
    MatrixFile_MAT4 h(file, rows, cols, name, type2prec<T>::prec);
    h.write_header();
    h.write_block(rows * cols, m.data());
  }
};

#endif  // FAST_PCA_FILE_MAT4_H_
