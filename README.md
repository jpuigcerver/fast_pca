fast_pca
========

A fast and memory efficient Principal Components Analysis software.

- Why fast?
- 
You can compute the mean, standard deviation and the eigenvalues and
eigenvectors of the covariance matrix with just a single pass through
your data. Plus, if your data is split into several matrices you
can do this work in parallel using a Mapreduce approach.

fast_pca uses BLAS + LAPACK to perform the required matrix operations efficiently. CMake 2.8 can detect several vendors implementing these interfaces like ATLAS, Goto, Intel MKL and many others.

- Why memory efficient?

You don't need to load the whole dataset into main memory not even a
single time. Plus, you can decide to use simple or double precision to
process your data.

- Why another PCA software?

There is plenty of software that can do the same job than this
one. However all tools that I knew forced you to load all your data
into main memory and then perform one or even multiple passes accross
your data, which was not feasible for some datasets that I was handling.
