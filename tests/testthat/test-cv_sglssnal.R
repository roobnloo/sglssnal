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
  expect_obj <- c(4.721, 4.721)
  names(expect_obj) <- c("primal objective", "dual objective")
  expect_equal(obj, expect_obj)
  expect_equal(result$info$iter, 8)
  expect_equal(result$info$nnz, 7)
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

test_that("lam_max scales with alphas", {
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

  result <- sglssnal::cv.sglssnal(
    A, ystar, grp, ind,
    alphas = 1, quietall = TRUE
  )

  lam_max <- result$cv_info$lambdas[1]

  result_01 <- sglssnal::cv.sglssnal(
    A, ystar, grp, ind,
    alphas = 0.1, quietall = TRUE
  )

  lam_max01 <- result_01$cv_info$lambdas[1]

  expect_equal(round(lam_max01, 2), round(lam_max * 10, 2))
})
