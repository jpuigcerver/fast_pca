fast_pca
========

A fast and memory efficient Principal Components Analysis software.

- **Why fast?**
  You can compute the mean, standard deviation and the eigenvalues and
  eigenvectors of the covariance matrix with just a single pass through
  your data. Plus, if your data is split into several matrices you
  can do this work in parallel using a Mapreduce approach.

  fast_pca uses BLAS + LAPACK to perform the required matrix operations
  efficiently. CMake 2.8 can detect several vendors implementing these
  interfaces like ATLAS, Goto, Intel MKL and many others.


- **Why memory efficient?**
  You don't need to load the whole dataset into main memory not even a
  single time. Plus, you can decide to use simple or double precision to
  process your data.


- **Why another PCA software?**
  There is plenty of software that can do the same job than this
  one. However all tools that I knew forced you to load all your data
  into main memory and then perform one or even multiple passes accross
  your data, which was not feasible for some datasets that I was handling.


### Minimum requirements
- C++ Compiler with capability for C++11. I tested with GCC 4.8.3 on
GNU/Linux and Clang 503.0.40 on Mac OS X.
- CMake 2.8
- LAPACK distribution. CMake will try detect automatically which LAPACK
distributions you have installed. If you have more than one, you can set
the CMake variable BLA_VENDOR to your favourite one. I tested the software
with the generic Netlib LAPACK implementation, ATLAS and
Mac OS X Accelerate framework.

### Installation
1. Download and unzip or clone the fast_pca repository into your
local machine.
2. Create a "build" directory in "fast_pca" directory:
```mkdir -p fast_pca/build && cd fast_pca/build```
3. Configure cmake:
```cmake ..```
4. Compile and install:
```make && make install```


### Basic usage:
- Compute mean, standard deviation and eigenvalues and eigenvectors of a
matrix:
```
fast_pca -C A.mat > pca.mat
```
This will generate a file containing 3 vectors (the mean M, the standard
deviation S and the eigenvalues D) and a matrix (the eigenvectors V). The
file is in the same format used by Octave.

The input data matrix can also be read from stdin.

- Project a matrix into a lower-dimensional space:
```
fast_pca -P -m pca.mat -q 2 A.mat > A.2d.mat
```
fast_pca reads the PCA data file and projects the data matrix A.mat into a
two-dimensional space. The resulting data is stored into A.2d.mat.

The input data matrix can also be read from stdin.

- Perform PCA and projection to preserve 95% of the variance using a
single call:
```
fast_pca -j 0.95 A.mat
```
fast_pca will perform the PCA and the projection on the same call.
Observe that the input data will be read twice: first to perform the PCA
and finally to perform de dimensionality reduction. The projected
data will be print on the standard output.


### Matrix formats:

For all formats, data samples are represented by row-vectors. Matrices are
read from the file and stored in memory assuming row-major order (all elements
from a given row are continuous).

I will assume from now on that your data matrix is an M x N matrix named A.

#### ASCII

This is the default format used by fast_pca and it consits only of the sequence
of rows which compose your data represented in ASCII characters. Elements of
the matrix are separated by whitespaces (one or many). Usually, one would
represent the different rows by different lines, but this is only optional.

Observe that since no header information is given regarding the number of
dimensions (columns) of your data, you will need to give that information
to the program through the -p option. By default, elements are read as
simple precision numbers, but you can interpret them as double (with -d).

```
A(1,1) A(1,2) ... A(1,N)
A(2,1) A(2,2) ... A(2,N)
...
A(M,1) A(M,2) ... A(M,N)
```

#### Binary

If you are working with very big datasets, ASCII representations may not be
convinient. The Binary format reads/writes all the matrix values in row-major
order into a file, with no header information.

Same as before, you will need to indicate the number of columns through the
-p option. Also, you also need to read the file using the same precision
that it was used to write it. Plus, no byte-ordering assumptions are made, thus
if you write a matrix in a machine working with Little-Endian ordering and then
try to use it on a Big-Endian machine, ugly things will happen.

#### Octave

fast_pca supports the ASCII Octave file format. The file representing our
previous MxN-matrix A, would be:

```
# name: A
# type: matrix
# rows: M
# columns: N
A(1,1) A(1,2) ... A(1,N)
A(2,1) A(2,2) ... A(2,N)
...
A(M,1) A(M,2) ... A(M,N)
```

The name of a matrix is an optional field (actually, fast_pca just ignores it).
Also, the header may be preceeded with other comment lines and they will be
ignored. However, note that the relative ordering between the "type", "rows"
and "columns" must be maintained, and the "columns" field must be the last
comment line before the matrix data.

#### VBosch

This format was used by other PRHLT tools and was added for backward
compatibility. The structure is very similar to the ASCII format, but adding
a header line indicating the number of rows and columns of the data matrix.

```
M N
A(1,1) A(1,2) ... A(1,N)
A(2,1) A(2,2) ... A(2,N)
...
A(M,1) A(M,2) ... A(M,N)
```