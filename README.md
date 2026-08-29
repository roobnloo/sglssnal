# Sparse-Group Lasso via Semismooth Newton Augmented Lagrangian (SGLSSNAL)

## Description

This R package solves the sparse-group lasso problem using second-order information via the Semismooth Newton Augmented Lagrangian method from [Zhang et al. (2020)](https://link.springer.com/article/10.1007/s10107-018-1329-6).

For a vector $x$ partitioned into $g$ non-overlapping groups and nonnegative weights $\{w_i\}$, define the penalty function

```math
\Phi(x) = \lambda_1\lVert x \rVert_1 + \lambda_2\sum_{i=1}^g w_i \lVert x_{(i)} \rVert_2.
```
The sparse-group lasso problem has the form
```math
\text{min}_x\ \frac{1}{2} \lVert Ax - b\rVert_2^2 + \Phi(x) \qquad \text{(P)}
```
while the dual problem is given by
```math
\begin{matrix}\max_{y,z}\; -\langle b, y \rangle - \frac{1}{2}\lVert y\rVert_2^2 - \Phi^\ast(z) \\ \mathrm{s.t.}\; A^\top y + z = 0 \qquad \text{(D)} \end{matrix}
```

where $\Phi^\ast$ denotes the convex conjugate.
Unlike first-order descent based methods which focus on $\mathrm{P}$, the SSNAL method uses second-order techniques to solve $\mathrm{D}$.

This package is designed for high-dimensional data with a (non-overlapping) group structure. The timing and accuracy of this method seems to be better than first-order descent methods in a large variety of problems.

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
A <- Matrix::rsparsematrix(n, p, density = 0.5, rand.x = rnorm)
ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

# group[j] is the group id of column j of A; here columns 1:20, 21:100,
# and 101:200 form three groups.
group <- rep(1:3, times = c(20, 80, 100))

result <- sglssnal::sglssnal(A, ystar, group, lambda = 2, alpha = 0.5)
print(result$info)
print(coef(result))
```
