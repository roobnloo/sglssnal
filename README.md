# Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian (SGLSSNAL)

## Description

This R package solves the sparse-group lasso problem using second-order information via the Semismooth Newton Augmented Lagrangian method from [Zhang et al. (2020)](https://link.springer.com/article/10.1007/s10107-018-1329-6).

For a vector $x$ split up into $g$ (possibly overlapping) groups and nonnegative weights $\{w_i\}$, define the penalty function 
$$\Phi(x) = \lambda_1\lVert x \rVert_1 + \lambda_2\sum_{i=1}^g w_i \lVert x_{(i)} \rVert_2.$$
The sparse-group problem has the form
$$\operatorname*{min}_x\; \frac{1}{2} \lVert Ax - b\rVert_2^2 + \Phi(x) \tag{P}$$
while the dual problem is given by
$$
\max_{y,z}\; -\langle b, y \rangle - \frac{1}{2}\lVert y\rVert_2^2 - \Phi^\ast(z)\\
\mathrm{s.t.}\; A^\top y + z = 0, \tag{D}
$$
where $\Phi^\ast$ denotes the convex conjugate.
Unlike first-order descent based methods which focus on $\mathrm{P}$, the SSNAL method uses second-order techniques to solve $\mathrm{D}$.

This package handles high-dimensional data with possibly overlapping group structures. The timing and accuracy of this method seems to be better than first-order descent methods in a large variety of problems.

The code was initially ported directly from [the original Matlab](https://github.com/YangjingZhang/SparseGroupLasso).

## Installation

```r
install.packages("devtools")
devtools::install_github("roobnloo/sglssnal")
```

## Example usage

```r
set.seed(11234)
n <- 100
p <- 200

bstar <- c(rnorm(20), rep(0, p - 20))
A <- matrix(rnorm(n * p), nrow = n)
A[sample(n * p, n * p / 2)] <- 0
A <- Matrix::Matrix(rnorm(n * p), nrow = n, sparse = TRUE)
ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

grp <- 1:p
ind <- matrix(c(1, 20, 1, 21, 100, 1, 101, 200, 1), nrow = 3)

result <- sglssnal::sglssnal(A, ystar, 0.1, 0.1, grp, ind)
print(result$info)
print(result$x)
```