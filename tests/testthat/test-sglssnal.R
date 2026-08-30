test_that("simple run of sglssnal", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1:3, times = c(20, 80, 100))

  result <- sglssnal::sglssnal(A, ystar, group, 2,
    alpha = 0.5, intercept = FALSE, standardize = FALSE
  )
  obj <- round(result$obj, 3)
  expect_obj <- matrix(c(19.358, 19.358), nrow = 2)
  rownames(expect_obj) <- c("primal objective", "dual objective")
  expect_equal(obj, expect_obj)
  expect_equal(result$info[[1]]$iter, 10)
  expect_equal(result$info[[1]]$nnz, 50)
})

test_that("sparse and dense run", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  A[sample(n * p, n * p / 2)] <- 0
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1:3, times = c(20, 80, 100))

  result_dense <- sglssnal::sglssnal(A, ystar, group, 2, alpha = 0.5)

  A_sp <- Matrix::Matrix(A, sparse = TRUE)
  result_sparse <- sglssnal::sglssnal(A_sp, ystar, group, 2, alpha = 0.5)

  obj <- sum(abs(round(result_dense$obj - result_sparse$obj, 10)))
  expect_equal(obj, 0)
})

test_that("intercept", {
  set.seed(231415)
  n <- 400
  p <- 2

  # Create data with a known intercept
  b0_true <- 5.0
  bstar <- rnorm(p, sd = 0.5)
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- b0_true + as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1, p)
  pfgroup <- 0

  result_intr <- sglssnal::sglssnal(
    A, ystar, group, 1,
    alpha = 0, pfgroup = pfgroup, stoptol = 1e-8, intercept = TRUE, standardize = FALSE
  )

  cm <- colMeans(A)
  Ac <- sweep(A, 2, cm, "-")
  yc <- ystar - mean(ystar)

  result_cent <- sglssnal::sglssnal(
    Ac, yc, group, 1,
    alpha = 0, pfgroup = pfgroup, stoptol = 1e-8, intercept = FALSE, standardize = FALSE
  )

  ints <- c(as.numeric(mean(ystar) - cm %*% result_cent$x), result_intr$x0)

  expect_true(
    abs(diff(ints)) < 0.01,
    "Intercept should be reasonably close to estimate based on centering"
  )
})

test_that("lambda path", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1:3, times = c(20, 80, 100))

  lambda <- seq(0.1, 2, length.out = 10)

  result <- sglssnal::sglssnal(A, ystar, group,
    lambda = lambda, alpha = 0.5, intercept = FALSE, standardize = FALSE
  )

  expect_equal(length(result$lambda), 10)
  expect_equal(result$lambda, sort(lambda, decreasing = TRUE))
  expect_equal(dim(result$x), c(p, 10))
  expect_equal(dim(result$y), c(n, 10))
  expect_equal(dim(result$z), c(p, 10))
})

test_that("group must match ncol(A)", {
  set.seed(231415)
  n <- 20
  p <- 10
  A <- matrix(rnorm(n * p), nrow = n)
  b <- rnorm(n)

  expect_error(
    sglssnal::sglssnal(A, b, group = 1:(p - 1)),
    "length\\(group\\) must be equal to ncol\\(A\\)"
  )
})

test_that("group must not contain NA", {
  set.seed(231415)
  n <- 20
  p <- 10
  A <- matrix(rnorm(n * p), nrow = n)
  b <- rnorm(n)
  group <- c(1:(p - 1), NA)

  expect_error(
    sglssnal::sglssnal(A, b, group = group),
    "group must not contain NA"
  )
})

test_that("pfgroup is indexed by sort(unique(group)), not first-appearance order", {
  set.seed(231415)
  n <- 50
  p <- 10
  A <- matrix(rnorm(n * p), nrow = n)
  bstar <- rnorm(p)
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  # Baseline: labels already in sorted == first-appearance order.
  group_a <- c(rep(1, 5), rep(2, 5))
  pfgroup_a <- c(1, 5)

  # Same column-to-weight assignment, but labels chosen so first-appearance
  # order (5, 1) differs from sort(unique(.)) order (1, 5). pfgroup is
  # supplied in sorted order: weight 5 for label 1 (cols 6:10), weight 1 for
  # label 5 (cols 1:5) -- matching group_a's cols 1:5 -> weight 1, cols
  # 6:10 -> weight 5.
  group_b <- c(rep(5, 5), rep(1, 5))
  pfgroup_b <- c(5, 1)

  result_a <- sglssnal::sglssnal(
    A, b, group_a,
    lambda = 1, alpha = 0.5, pfgroup = pfgroup_a,
    intercept = FALSE, standardize = FALSE
  )
  result_b <- sglssnal::sglssnal(
    A, b, group_b,
    lambda = 1, alpha = 0.5, pfgroup = pfgroup_b,
    intercept = FALSE, standardize = FALSE
  )

  expect_equal(as.matrix(result_a$x), as.matrix(result_b$x))
})

test_that("standardize = TRUE matches manually standardizing and unscaling", {
  set.seed(1)
  n <- 100
  p <- 20

  A <- matrix(rnorm(n * p), n, p)
  bstar <- c(rnorm(5), rep(0, p - 5))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:4, each = 5)

  fit_std <- sglssnal(A, b, group,
    lambda = 1, alpha = 0.5, intercept = FALSE, standardize = TRUE
  )

  # standardize = TRUE fits on A scaled to unit column norm, then unscales
  # the coefficients back to A's original units -- replicate that manually
  # with standardize = FALSE and check they agree.
  A_sd <- sqrt(colSums(A^2))
  A_scaled <- sweep(A, 2, A_sd, "/")
  fit_manual <- sglssnal(A_scaled, b, group,
    lambda = 1, alpha = 0.5, intercept = FALSE, standardize = FALSE
  )

  expect_equal(as.numeric(fit_std$x), as.numeric(fit_manual$x) / A_sd, tolerance = 1e-6)
})

test_that("pfgroup weights shrink a more heavily penalized group harder", {
  set.seed(3)
  n <- 100
  p <- 10

  A <- matrix(rnorm(n * p), n, p)
  bstar <- rep(1, p)
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:2, each = 5)

  # Both groups carry identical true signal, so with equal weights their
  # fitted group norms should be nearly equal; penalizing group 2 harder
  # should shrink it well below group 1's norm.
  fit_equal <- sglssnal(A, b, group,
    lambda = 5, alpha = 0, pfgroup = c(1, 1),
    intercept = FALSE, standardize = FALSE
  )
  fit_weighted <- sglssnal(A, b, group,
    lambda = 5, alpha = 0, pfgroup = c(1, 8),
    intercept = FALSE, standardize = FALSE
  )

  norm1_equal <- sqrt(sum(fit_equal$x[1:5, 1]^2))
  norm2_equal <- sqrt(sum(fit_equal$x[6:10, 1]^2))
  norm1_weighted <- sqrt(sum(fit_weighted$x[1:5, 1]^2))
  norm2_weighted <- sqrt(sum(fit_weighted$x[6:10, 1]^2))

  expect_equal(norm1_equal, norm2_equal, tolerance = 0.05)
  expect_lt(norm2_weighted, norm2_equal)
  expect_gt(norm1_weighted, norm2_weighted)
})

test_that("alpha = 1 (pure lasso) matches the closed-form soft-threshold solution", {
  # With alpha = 1, lambda2 = 0, so the group penalty vanishes and this
  # reduces to plain lasso regardless of grouping. For an orthonormal
  # design (A'A = I), the lasso solution has the closed form
  # soft_threshold(A'b, lambda1) -- a strong, tuning-free correctness check.
  set.seed(7)
  n <- 50
  p <- 10

  A <- qr.Q(qr(matrix(rnorm(n * p), n, p)))
  bstar <- c(2, -1.5, rep(0, p - 2))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.05))
  group <- rep(1:5, each = 2)
  lambda_val <- 0.3

  fit <- sglssnal(A, b, group,
    lambda = lambda_val, alpha = 1,
    intercept = FALSE, standardize = FALSE, stoptol = 1e-10
  )

  z <- as.numeric(crossprod(A, b))
  expected <- sign(z) * pmax(abs(z) - lambda_val, 0)

  expect_equal(as.numeric(fit$x), expected, tolerance = 1e-6)
})

test_that("stopopt variants 1, 3, 4 agree with the default (2)", {
  set.seed(11)
  n <- 80
  p <- 15

  A <- matrix(rnorm(n * p), n, p)
  bstar <- c(rnorm(5), rep(0, p - 5))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:3, each = 5)

  base <- sglssnal(A, b, group,
    lambda = 0.5, alpha = 0.5, stopopt = 2L,
    intercept = FALSE, standardize = FALSE, stoptol = 1e-8
  )

  for (so in c(1L, 3L, 4L)) {
    fit <- sglssnal(A, b, group,
      lambda = 0.5, alpha = 0.5, stopopt = so,
      intercept = FALSE, standardize = FALSE, stoptol = 1e-8
    )
    expect_equal(as.numeric(fit$x), as.numeric(base$x), tolerance = 1e-6)
  }
})

test_that("user-supplied Lip/y0/z0/x0 warm-start an already-converged solution", {
  set.seed(21)
  n <- 60
  p <- 10

  A <- matrix(rnorm(n * p), n, p)
  bstar <- c(rnorm(3), rep(0, p - 3))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:2, each = 5)

  fit0 <- sglssnal(A, b, group,
    lambda = 0.5, alpha = 0.5, intercept = FALSE, standardize = FALSE
  )
  fit_warm <- sglssnal(A, b, group,
    lambda = 0.5, alpha = 0.5, intercept = FALSE, standardize = FALSE,
    Lip = sglssnal:::compute_lip(A),
    y0 = fit0$y[, 1], z0 = fit0$z[, 1], x0 = as.numeric(fit0$x[, 1])
  )

  # `iter` counts the augmented-Lagrangian outer loop (sglssnal_main), not
  # the semismooth-Newton inner iterations (tracked separately as
  # `itersub`). Starting from an already-converged iterate should take a
  # single outer iteration to re-confirm convergence.
  expect_equal(fit_warm$info[[1]]$iter, 1)
  expect_equal(as.numeric(fit_warm$x), as.numeric(fit0$x), tolerance = 1e-6)
})

test_that("mismatched y0/z0/x0 dimensions raise a clear error", {
  set.seed(21)
  n <- 60
  p <- 10
  A <- matrix(rnorm(n * p), n, p)
  b <- rnorm(n)
  group <- rep(1:2, each = 5)

  expect_error(
    sglssnal(A, b, group, y0 = rep(0, n), z0 = rep(0, p), x0 = rep(0, p - 1)),
    "y0, z0, and x0 must have dimensions"
  )
})

test_that("invalid arguments are rejected with informative errors", {
  set.seed(21)
  n <- 60
  p <- 10
  A <- matrix(rnorm(n * p), n, p)
  b <- rnorm(n)
  group <- rep(1:2, each = 5)

  expect_error(sglssnal(A, b, group, verbose = 9L), "verbose must be one of 0, 1, or 2")
  expect_error(sglssnal(A, b, group, stopopt = 5L), "stopopt must be one of 1, 2, 3, or 4")
  expect_error(sglssnal(A, b[-1], group), "nrow\\(A\\) must be equal to length\\(b\\)")
  expect_error(sglssnal(A, b, group, alpha = 1.5), "alpha must be in \\[0, 1\\]")
  expect_error(sglssnal(A, b, group, alpha = c(0.1, 0.2)), "length\\(alpha\\) must be 1")
  expect_error(sglssnal(A, b, group, lambda = -1), "lambda values must be positive")
  expect_error(sglssnal(A, b, group, maxit = 0), "maxit must be a positive integer")
  expect_error(sglssnal(A, b, group, maxit = 3.5), "maxit must be a positive integer")
  expect_error(sglssnal(A, b, group, nlambda = 3.5), "nlambda must be a positive integer")
  expect_error(sglssnal(A, b, group, maxit = c(5, 10)), "length\\(maxit\\) must be 1")
  expect_error(sglssnal(A, b, group, nlambda = c(5, 10)), "nlambda must be a positive integer")
  expect_error(sglssnal(A, b, group, stoptol = -1), "stoptol must be a positive number")
  expect_error(
    sglssnal(A, b, group, pfgroup = c(1, 1, 1)),
    "length\\(pfgroup\\) must equal the number of unique groups"
  )

  A_empty <- matrix(nrow = n, ncol = 0)
  expect_error(
    sglssnal(A_empty, b, integer(0), lambda = 1),
    "A must have at least one row and one column"
  )

  A_na <- A
  A_na[1, 1] <- NA
  expect_error(
    sglssnal(A_na, b, group, lambda = 1),
    "A must not contain missing, NaN, or infinite values"
  )
  A_inf <- A
  A_inf[1, 1] <- Inf
  expect_error(
    sglssnal(A_inf, b, group, lambda = 1),
    "A must not contain missing, NaN, or infinite values"
  )

  b_na <- b
  b_na[1] <- NA
  expect_error(
    sglssnal(A, b_na, group, lambda = 1),
    "b must not contain missing, NaN, or infinite values"
  )

  A_sp_inf <- Matrix::Matrix(A, sparse = TRUE)
  A_sp_inf[1, 1] <- Inf
  expect_error(
    sglssnal(A_sp_inf, b, group, lambda = 1),
    "A must not contain missing, NaN, or infinite values"
  )
})

test_that("large active-set problem exercises the pcg linear solver path", {
  # n = 150 observations, p = 9000 predictors, alpha = 0 (pure group lasso,
  # no l1 thresholding) with a small lambda relative to signal keeps most
  # coordinates active, driving the internal active-set size (`nnz`, what
  # conjgrad_linsolver's dispatcher calls `density`) well past the 8000
  # threshold that routes n=150 to solver 3 (pcg) instead of the dense
  # direct solve -- see test-conjgrad-linsolver.R for the actual numeric
  # correctness check of that path; this is an end-to-end smoke test that
  # it's wired up and doesn't crash/hang through the public API.
  set.seed(42)
  n <- 150
  p <- 9000
  A <- matrix(rnorm(n * p), nrow = n)
  bstar <- rnorm(p, sd = 0.5)
  y <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:10, each = p / 10)

  result <- sglssnal::sglssnal(
    A, y, group,
    lambda = 0.01, alpha = 0, stoptol = 1e-2, maxit = 50, intercept = FALSE, standardize = FALSE
  )

  expect_true(all(is.finite(result$x@x)))
  expect_equal(result$info[[1]]$msg, "converged")
  expect_gt(result$info[[1]]$nnz, 8000)
})

test_that("auto lambda path is unaffected by a large constant offset in b", {
  # With intercept = TRUE, adding a constant to b must not change the fit
  # at all -- the intercept absorbs it entirely.
  set.seed(1)
  n <- 100
  p <- 20
  A <- matrix(rnorm(n * p), n, p)
  bstar <- c(rnorm(5), rep(0, p - 5))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:4, each = 5)

  fit1 <- sglssnal(A, b, group, nlambda = 10, alpha = 0.5)
  fit2 <- sglssnal(A, b + 1000, group, nlambda = 10, alpha = 0.5)

  expect_equal(fit1$lambda, fit2$lambda, tolerance = 1e-6)
  expect_equal(as.matrix(fit1$x), as.matrix(fit2$x), tolerance = 1e-6)
  expect_false(all(as.matrix(fit1$x) == 0))
})

test_that("auto lambda path works for sparse A without Matrix attached", {
  # sglssnal() must work for sparse A even when Matrix is only imported,
  # not attached via library(Matrix).
  skip_if("Matrix" %in% .packages(), "Matrix must not be attached for this check")

  set.seed(1)
  n <- 50
  p <- 20
  A <- Matrix::rsparsematrix(n, p, density = 0.3, rand.x = rnorm)
  bstar <- c(rnorm(5), rep(0, p - 5))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:4, each = 5)

  fit <- sglssnal(A, b, group, nlambda = 10, alpha = 0.5)
  expect_equal(dim(fit$x), c(p, 10))
})
