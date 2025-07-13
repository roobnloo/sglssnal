test_that("run with custom lambda", {
  set.seed(231415)
  n <- 100
  p <- 20

  bstar <- c(rnorm(5), rep(0, p - 5))
  A <- matrix(rnorm(n * p), nrow = n)
  A[sample(n * p, n * p / 2)] <- 0
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  grp <- 1:p
  ind <- matrix(c(1, 10, 11, 15, 16, 20), nrow = 2)

  lambda <- seq(3, 0.01, length.out = 10)

  result <- sglssnal::cv.sglssnal(
    A, ystar, grp, ind,
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

  grp <- 1:p
  ind <- matrix(c(1, 10, 11, 15, 16, 20), nrow = 2)

  nlambda <- 10
  result <- sglssnal::cv.sglssnal(
    A, ystar, grp, ind,
    nlambda = nlambda, alpha = 0.5,
    quietall = TRUE, intercept = FALSE, standardize = FALSE
  )

  single <- sglssnal::sglssnal(
    A, ystar, grp, ind,
    nlambda = nlambda, alpha = 0.5,
    intercept = FALSE, standardize = FALSE, printyes = FALSE
  )

  expect_true(all.equal(result$x, single$x))
})
