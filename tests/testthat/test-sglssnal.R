test_that("simple run of sglssnal", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1:3, times = c(20, 80, 100))

  result <- sglssnal::sglssnal(A, ystar, group, 2,
    alpha = 0.5, printmain = FALSE, intercept = FALSE, standardize = FALSE
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

  result_dense <- sglssnal::sglssnal(A, ystar, group, 2, alpha = 0.5, printmain = FALSE)

  A_sp <- Matrix::Matrix(A, sparse = TRUE)
  result_sparse <- sglssnal::sglssnal(A_sp, ystar, group, 2, alpha = 0.5, printmain = FALSE)

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
    alpha = 0, pfgroup = pfgroup, stoptol = 1e-8,
    printmain = FALSE, intercept = TRUE, standardize = FALSE
  )

  cm <- colMeans(A)
  Ac <- sweep(A, 2, cm, "-")
  yc <- ystar - mean(ystar)

  result_cent <- sglssnal::sglssnal(
    Ac, yc, group, 1,
    alpha = 0, pfgroup = pfgroup, stoptol = 1e-8,
    printmain = FALSE, intercept = FALSE, standardize = FALSE
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
    lambda = lambda, alpha = 0.5,
    printmain = FALSE, intercept = FALSE, standardize = FALSE
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
    sglssnal::sglssnal(A, b, group = 1:(p - 1), printmain = FALSE),
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
    sglssnal::sglssnal(A, b, group = group, printmain = FALSE),
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
    intercept = FALSE, standardize = FALSE, printmain = FALSE
  )
  result_b <- sglssnal::sglssnal(
    A, b, group_b,
    lambda = 1, alpha = 0.5, pfgroup = pfgroup_b,
    intercept = FALSE, standardize = FALSE, printmain = FALSE
  )

  expect_equal(as.matrix(result_a$x), as.matrix(result_b$x))
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
    lambda = 0.01, alpha = 0, stoptol = 1e-2, maxit = 50,
    printmain = FALSE, intercept = FALSE, standardize = FALSE
  )

  expect_true(all(is.finite(result$x@x)))
  expect_equal(result$info[[1]]$msg, "converged")
  expect_gt(result$info[[1]]$nnz, 8000)
})
