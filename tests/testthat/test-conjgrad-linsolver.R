test_that("pcg (solver 3) matches direct dense solve (solver 1) - sparse A", {
  set.seed(1)
  n <- 150
  p <- 40
  A <- matrix(rnorm(n * p), nrow = n)
  A_sp <- Matrix::Matrix(A, sparse = TRUE)
  rhs <- rnorm(n)
  u <- rnorm(p)
  group <- rep(1:4, each = 10)
  gs_list <- sglssnal:::group_structure(p, group, rep(1, 4))
  lam1 <- 0.1
  lam2 <- 0.2
  sig <- 1.0
  par <- list(tol = 1e-10, maxit = 5000, precond = 0)

  # density = 0L keeps n <= 1000 routed to solver 1 (dense direct, ground
  # truth); density = 9000L forces solver 3 via the (n > 100 && density >
  # 8000) branch regardless of the true active-set size.
  r1 <- sglssnal:::conjgrad_linsolver_interface(
    A_sp, rhs, u, lam1, lam2, gs_list, density = 0L, sig = sig, par = par
  )
  r3 <- sglssnal:::conjgrad_linsolver_interface(
    A_sp, rhs, u, lam1, lam2, gs_list, density = 9000L, sig = sig, par = par
  )

  expect_equal(as.numeric(r1$dy), as.numeric(r3$dy), tolerance = 1e-8)
  expect_equal(r1$solve_ok, 1L)
  expect_equal(r3$solve_ok, 1L)
  expect_true(length(r3$resnrm) > 1) # confirms psqmr actually iterated
})

test_that("pcg (solver 3) matches direct dense solve (solver 1) - dense A", {
  set.seed(2)
  n <- 150
  p <- 40
  A <- matrix(rnorm(n * p), nrow = n)
  rhs <- rnorm(n)
  u <- rnorm(p)
  group <- rep(1:4, each = 10)
  gs_list <- sglssnal:::group_structure(p, group, rep(1, 4))
  lam1 <- 0.1
  lam2 <- 0.2
  sig <- 1.0
  par <- list(tol = 1e-10, maxit = 5000, precond = 0)

  r1 <- sglssnal:::conjgrad_linsolver_interface_dense(
    A, rhs, u, lam1, lam2, gs_list, density = 0L, sig = sig, par = par
  )
  r3 <- sglssnal:::conjgrad_linsolver_interface_dense(
    A, rhs, u, lam1, lam2, gs_list, density = 9000L, sig = sig, par = par
  )

  expect_equal(as.numeric(r1$dy), as.numeric(r3$dy), tolerance = 1e-8)
  expect_true(length(r3$resnrm) > 1)
})

test_that("solver 3 handles the no-active-groups (identity) case", {
  set.seed(3)
  n <- 150
  p <- 40
  A <- Matrix::Matrix(matrix(rnorm(n * p), nrow = n), sparse = TRUE)
  rhs <- rnorm(n)
  u <- rep(0.01, p) # small enough that proximal_l1(u, lam1) is all zero
  group <- rep(1:4, each = 10)
  gs_list <- sglssnal:::group_structure(p, group, rep(1, 4))
  lam1 <- 0.1
  lam2 <- 0.2
  sig <- 1.0
  par <- list(tol = 1e-10, maxit = 5000, precond = 0)

  r3 <- sglssnal:::conjgrad_linsolver_interface(
    A, rhs, u, lam1, lam2, gs_list, density = 9000L, sig = sig, par = par
  )

  expect_equal(as.numeric(r3$dy), rhs)
  expect_equal(r3$solve_ok, 1L)
  expect_equal(length(r3$resnrm), 1L)
})
