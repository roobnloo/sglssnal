test_that("simple run of sglssnal", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  grp <- 1:p
  ind <- matrix(c(1, 20, 21, 100, 101, 200), nrow = 2)

  result <- sglssnal::sglssnal(A, ystar, grp, ind, 2,
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

  grp <- 1:p
  ind <- matrix(c(1, 20, 21, 100, 101, 200), nrow = 2)

  result_dense <- sglssnal::sglssnal(A, ystar, grp, ind, 2, alpha = 0.5, printmain = FALSE)

  A_sp <- Matrix::Matrix(A, sparse = TRUE)
  result_sparse <- sglssnal::sglssnal(A_sp, ystar, grp, ind, 2, alpha = 0.5, printmain = FALSE)

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

  grp <- 1:p
  ind <- matrix(c(1, p), nrow = 2)
  pfgroup <- 0

  result_intr <- sglssnal::sglssnal(
    A, ystar, grp, ind, 1,
    alpha = 0, pfgroup = pfgroup, stoptol = 1e-8,
    printmain = FALSE, intercept = TRUE, standardize = FALSE
  )

  cm <- colMeans(A)
  Ac <- sweep(A, 2, cm, "-")
  yc <- ystar - mean(ystar)

  result_cent <- sglssnal::sglssnal(
    Ac, yc, grp, ind, 1,
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

  grp <- 1:p
  ind <- matrix(c(1, 20, 21, 100, 101, 200), nrow = 2)

  lambda <- seq(0.1, 2, length.out = 10)

  result <- sglssnal::sglssnal(A, ystar, grp, ind,
    lambda = lambda, alpha = 0.5,
    printmain = FALSE, intercept = FALSE, standardize = FALSE
  )

  expect_equal(length(result$lambda), 10)
  expect_equal(result$lambda, sort(lambda, decreasing = TRUE))
  expect_equal(dim(result$x), c(p, 10))
  expect_equal(dim(result$y), c(n, 10))
  expect_equal(dim(result$z), c(p, 10))
})
