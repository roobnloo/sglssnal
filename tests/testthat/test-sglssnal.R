test_that("simple run of sglssnal", {
  set.seed(231415)
  n <- 100
  p <- 200

  bstar <- c(rnorm(20), rep(0, p - 20))
  A <- matrix(rnorm(n * p), nrow = n)
  A[sample(n * p, n * p / 2)] <- 0
  A <- matrix(rnorm(n * p), nrow = n)
  ystar <- as.numeric(A %*% bstar + rnorm(n, sd = 0.1))

  grp <- 1:p
  ind <- matrix(c(1, 20, 21, 100, 101, 200), nrow = 2)

  result <- sglssnal::sglssnal(A, ystar, grp, ind, 2, 0.5, printyes = FALSE)
  obj <- round(result$obj, 3)
  expect_obj <- c(19.312, 19.312)
  names(expect_obj) <- c("primal objective", "dual objective")
  expect_equal(obj, expect_obj)
  expect_equal(result$info$iter, 10)
  expect_equal(result$info$nnz, 47)
})
