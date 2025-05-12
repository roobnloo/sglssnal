test_that("simple run of sglssnal", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  grp <- 1:p
  ind <- matrix(c(1, 20, 21, 100, 101, 200), nrow = 2)

  result <- sglssnal::sglssnal(A, ystar, grp, ind, 2, 0.5,
    printyes = FALSE,
    intercept = FALSE, standardize = FALSE
  )
  obj <- round(result$obj, 3)
  expect_obj <- c(19.358, 19.358)
  names(expect_obj) <- c("primal objective", "dual objective")
  expect_equal(obj, expect_obj)
  expect_equal(result$info$iter, 10)
  expect_equal(result$info$nnz, 50)
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

  result_dense <- sglssnal::sglssnal(A, ystar, grp, ind, 2, 0.5, printyes = FALSE)

  A_sp <- Matrix::Matrix(A, sparse = TRUE)
  result_sparse <- sglssnal::sglssnal(A_sp, ystar, grp, ind, 2, 0.5, printyes = FALSE)

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
    A, ystar, grp, ind, 1, 0,
    pfgroup = pfgroup, stoptol = 1e-8,
    printyes = FALSE, intercept = TRUE, standardize = FALSE
  )

  cm <- colMeans(A)
  Ac <- sweep(A, 2, cm, "-")
  yc <- ystar - mean(ystar)

  result_cent <- sglssnal::sglssnal(
    Ac, yc, grp, ind, 1, 0,
    pfgroup = pfgroup, stoptol = 1e-8,
    printyes = FALSE, intercept = FALSE, standardize = FALSE
  )

  ints <- c(mean(ystar) - cm %*% result_cent$x, result_intr$x0)

  expect_true(
    abs(diff(ints)) < 0.01,
    "Intercept should be reasonably close to estimate based on centering"
  )
})
