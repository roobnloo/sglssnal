test_that("simple cv run", {
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

  result <- sglssnal::cv.sglssnal(A, ystar, grp, ind, quietall = TRUE)
  obj <- round(result$obj, 3)
  expect_obj <- c(8.9090, 8.9090)
  names(expect_obj) <- c("primal objective", "dual objective")
  expect_equal(obj, expect_obj)
  expect_equal(result$info$iter, 7)
  expect_equal(result$info$nnz, 5)
})

test_that("run with custom grid", {
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

  lambdas <- c(1, 0.5, 0.25, 0.1)
  alphas <- c(0.1, 0.5, 0.9)

  result <- sglssnal::cv.sglssnal(
    A, ystar, grp, ind,
    alphas = c(0.1, 0.5, 0.9), lambdas = c(1, 0.5, 0.25, 0.1), quietall = TRUE
  )
  obj <- round(result$obj, 3)
  expect_equal(obj, c("primal objective" = 5.977, "dual objective" = 5.977))
  expect_equal(result$info$iter, 9)
  expect_equal(result$info$nnz, 6)
  expect_equal(result$cv_info$lambdas, lambdas)
  expect_equal(result$cv_info$alphas, alphas)
  expect_equal(result$cv_info$cv_idx, c(cv_lambda_id = 1, cv_alpha_id = 3))
})
