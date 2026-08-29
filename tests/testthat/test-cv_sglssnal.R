test_that("run with custom lambda", {
  set.seed(231415)
  n <- 100
  p <- 20

  bstar <- c(rnorm(5), rep(0, p - 5))
  A <- matrix(rnorm(n * p), nrow = n)
  A[sample(n * p, n * p / 2)] <- 0
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1:3, times = c(10, 5, 5))

  lambda <- seq(3, 0.01, length.out = 10)

  result <- sglssnal::cv.sglssnal(
    A, ystar, group,
    lambda = lambda, alpha = 0.5, intercept = FALSE, standardize = FALSE
  )
  obj <- round(result$obj, 3)
  expect_equal(dim(obj), c(2, length(lambda)))
  expect_equal(result$cv_info$lambda, lambda)
  expect_equal(result$cv_info$alpha, 0.5)
  expect_equal(length(result$cv_info$cvm), length(lambda))
  expect_equal(result$cv_info$cv_lambda_id, 7)
})

test_that("run with default lambda", {
  set.seed(231415)
  n <- 100
  p <- 20

  bstar <- c(rnorm(5), rep(0, p - 5))
  A <- matrix(rnorm(n * p), nrow = n)
  A[sample(n * p, n * p / 2)] <- 0
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  group <- rep(1:3, times = c(10, 5, 5))

  nlambda <- 10
  result <- sglssnal::cv.sglssnal(
    A, ystar, group,
    nlambda = nlambda, alpha = 0.5, intercept = FALSE, standardize = FALSE
  )

  single <- sglssnal::sglssnal(
    A, ystar, group,
    nlambda = nlambda, alpha = 0.5,
    intercept = FALSE, standardize = FALSE
  )

  expect_true(all.equal(result$x, single$x))
})

test_that("custom foldid runs and produces valid cv_info", {
  set.seed(5)
  n <- 60
  p <- 10

  A <- matrix(rnorm(n * p), n, p)
  bstar <- c(rnorm(3), rep(0, p - 3))
  b <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))
  group <- rep(1:2, each = 5)

  foldid <- as.integer(rep(1:3, length.out = n))
  cvfit <- cv.sglssnal(A, b, group,
    nlambda = 5, foldid = foldid,
    intercept = FALSE, standardize = FALSE
  )

  expect_equal(length(cvfit$cv_info$cvm), 5)
  expect_true(cvfit$cv_info$cv_lambda_id %in% seq_len(5))
  expect_true(all(is.finite(cvfit$cv_info$cvm)))
})

test_that("invalid arguments are rejected with informative errors", {
  set.seed(5)
  n <- 60
  p <- 10
  A <- matrix(rnorm(n * p), n, p)
  b <- rnorm(n)
  group <- rep(1:2, each = 5)

  expect_error(cv.sglssnal(A, b, group, verbose = 9L), "verbose must be one of 0, 1, or 2")
  expect_error(cv.sglssnal(A, b, group, alpha = 1.5), "alpha must be in \\[0, 1\\]")
  expect_error(cv.sglssnal(A, b, group, alpha = c(0.1, 0.2)), "length\\(alpha\\) must be 1")
  expect_error(
    cv.sglssnal(A, b, group, foldid = as.integer(c(1, 1, 1))),
    "length\\(foldid\\) must be equal to nrow\\(A\\)"
  )
})
