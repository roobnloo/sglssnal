# Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian (SGLSSNAL)

## Description

This R package solves the sparse-group lasso problem using second-order information via the Semismooth Newton Augmented Lagrangian method from [Zhang et al. (2020)](https://link.springer.com/article/10.1007/s10107-018-1329-6).
It is designed to handle high-dimensional data with group structures, providing efficient and accurate solutions.
The timing and accuracy seems to be better than first-order descent methods in a large variety of problems.

The code was initially ported directly from [the original Matlab](https://github.com/YangjingZhang/SparseGroupLasso).

## Installation

```r
install.packages("devtools")
devtools::install_github("roobnloo/sglssnal")
```