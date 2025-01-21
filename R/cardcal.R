# Computes number of non-zero elements of x up to a tolerance of r.

cardcal <- function(x, r = 0.999) {
  n <- length(x)
  normx1 <- sum(abs(x))
  idx <- order(abs(x), decreasing = TRUE)
  absx <- abs(x)[idx]

  cumsum_absx <- cumsum(absx)
  k <- which(cumsum_absx >= r * normx1)[1]

  xnew <- rep(0, n)
  idxnew <- idx[1:k]
  xnew[idxnew] <- x[idxnew]

  return(list(k = k, xnew = xnew))
}
