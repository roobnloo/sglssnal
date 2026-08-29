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
    lambda = lambda, alpha = 0.5,
    quietall = TRUE, intercept = FALSE, standardize = FALSE
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
    nlambda = nlambda, alpha = 0.5,
    quietall = TRUE, intercept = FALSE, standardize = FALSE
  )

  single <- sglssnal::sglssnal(
    A, ystar, group,
    nlambda = nlambda, alpha = 0.5,
    intercept = FALSE, standardize = FALSE, printmain = FALSE
  )

  expect_true(all.equal(result$x, single$x))
})
