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
2. Create a "build" directory in fast_pca: 
```mkdir -p fast_pca/build && cd fast_pca_build```
3. Configure cmake: 
```cmake ..```
4. Compile and install:
```make && make_install```


### Basic usage:
- Compute mean, standard deviation and eigenvalues and eigenvectors of a
matrix:
```
fast_pca -C -e A.eigvec.mat -g A.eigval.mat -m A.mean.mat -s A.stddev.mat A.mat
```
This will generate the matrices A.eigvec.mat, A.eigval.mat (eigenvectors
and eigenvalues of the zero-mean data), A.mean.mat (mean of each 
dimension) and A.stddev.mat (standard deviation of each dimension).

The input data matrix can also be read from stdin.

- Project a matrix into a lower-dimensional space:
```
fast_pca -P -e A.eigvec.mat -g A.eigval.mat -m A.mean.mat -s A.stddev.mat -q 2 A.mat > A.2d.mat
``` 
fast_pca reads the eigenvectors, eigenvalues, mean and standard deviation
matrices and projects the data matrix A.mat into a two-dimensional space.
The resulting data is stored into A.2d.mat.

The input data matrix can also be read from stdin.

- Perform PCA and projection to preserve 95% of the variance using a
single call:
```
fast_pca -j 0.95 A.mat > A.2d.mat
```
fast_pca will perform the PCA and the projection on the same call.
Observe that the input data will be read twice: first to perform the PCA
and finally to perform de dimensionality reduction.


### Matrix formats:
- Data samples are represented by row-vectors. Matrices are read and 
stored in row-major order.
- For an M x N matrix A, the default format (called "text") is an ASCII 
file with the following format:
```
M N
A(1,1) A(1,2) ... A(1,N)
A(2,1) A(2,2) ... A(2,N)
...
A(M,1) A(M,2) ... A(M,N)
```
- Two additional formats which omit the header information are supported: 
"ascii" and "binary". In the binary version, the raw data matrix is stored. 
Observe that, when reading such matrices, you will need to specify the 
number of dimensions (columns) and the precision (simple vs. double) of 
the numeric data. These formats are convinient when reading streams of 
data, where the number of rows is not know in advance.
