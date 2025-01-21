conjgrad_linsolver <- function(A, rhs, u, c, P, par) {
  n <- length(rhs)
  solver <- 1 # 1: direct, 2: direct woodbury-formula, 3: pcg
  density <- par$nnz # l_1 sparsity
  dn <- 10000
  sig <- par$sigma

  if (n <= 1000) {
    solver <- 1
  }

  if (density <= n && density <= dn) {
    solver <- 2
  }

  # if ((n > 5000 && density >= 1000) || (n > 2000 && density > 5000) || (n > 100 && density > 8000)) {
  # TODO implement this heavy duty solver
  # solver <- 3
  # }

  resnrm <- NULL
  if (solver == 1) {
    mat_ssn_out <- mat_ssn(u, A, c, P, sig)
    # if (n <= 1000) {
    dy <- solve(mat_ssn_out$V, rhs)
    # } else {
    #   LAAt <- mycholAAt(mat_ssn_out$V, n)
    #   dy <- mylinsysolve(LAAt, rhs)
    # }
    resnrm <- 0
    solve_ok <- 1
  }

  if (solver == 2) {
    result <- mat2_ssn(u, A, c, P, sig)
    V2 <- result$V2
    D <- result$D
    id_yes <- result$id_yes

    if (id_yes) {
      dy <- rhs
    } else {
      rhstmp <- as.matrix(crossprod(D, rhs))
      dy <- solve(V2, rhstmp)
      dy <- D %*% dy
      dy <- rhs - dy
    }
    resnrm <- 0
    solve_ok <- 1
  }

  # if (solver == 3) {
  #   if (Ayes) {
  #     result <- matvecD(u, Ainput$A, c, P, sig)
  #     D <- result$D
  #     id_yes <- result$id_yes

  #     if (id_yes) {
  #       dy <- rhs
  #     } else {
  #       rhstmp <- crossprod(D, rhs)
  #       Afun <- function(x) x + crossprod(D, D %*% x)
  #       result <- psqmrGL(Afun, rhstmp, par)
  #       dy <- result$dy
  #       resnrm <- result$resnrm
  #       solve_ok <- result$solve_ok
  #       dy <- D %*% dy
  #       dy <- rhs - dy
  #     }
  #   } else {
  #     Afun <- function(x) matvecA(u, Ainput$Amap, Ainput$ATmap, c, P, sig, x)
  #     result <- psqmrGL(Afun, rhs, par)
  #     dy <- result$dy
  #     resnrm <- result$resnrm
  #     solve_ok <- result$solve_ok
  #   }
  # }

  return(list(dy = dy, resnrm = resnrm, solve_ok = solve_ok))
}
