fast_pca
========

A fast and memory efficient Principal Components Analysis software.

## Motivation

If you are familiar with Principal Components Analysis (PCA) you would
know that there are plenty of existing alternatives to reduce the
dimensionality of your favorite dataset. You can even use standard
numerical computing tools like Numpy + Scipy, MATLAB, etc.

The standard way of computing the PCA of a matrix requires that you
load the matrix into main memory and then compute the covariance
matrix and its eigenvectors.

However, sometimes it is not feasible to load the whole matrix into
main memory. In those cases, the previous software may be useless. 

Then, fast_pca may be a good solution for you.

## Usage

