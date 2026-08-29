test_that("compute_lip agrees between dense and sparse representations", {
  set.seed(1)
  n <- 50
  p <- 30
  A_dense <- matrix(rnorm(n * p), n, p)
  A_sparse <- Matrix::Matrix(A_dense, sparse = TRUE)

  lip_dense <- sglssnal:::compute_lip(A_dense, n)
  lip_sparse <- sglssnal:::compute_lip(A_sparse, n)

  expect_equal(lip_dense, lip_sparse, tolerance = 1e-8)
})

test_that("compute_lip matches a direct eigen() computation", {
  set.seed(2)
  n <- 20
  p <- 15
  A <- matrix(rnorm(n * p), n, p)

  lip <- sglssnal:::compute_lip(A, n)
  expected <- max(eigen(tcrossprod(A), symmetric = TRUE, only.values = TRUE)$values)

  expect_equal(as.numeric(lip), expected, tolerance = 1e-8)
})
