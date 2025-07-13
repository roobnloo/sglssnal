#' Extract estimated coefficients from a `sglssnal` object
#' @param object Fitted object of class `sglssnal`.
#' @param ... Additional arguments passed to or from other methods.
#' @return A numeric matrix of estimated coefficients.
#' @method coef sglssnal
#' @export
coef.sglssnal <- function(object, ...) {
  x0 <- matrix(object$x0, nrow = 1)
  rownames(x0) <- "(Intercept)"
  x <- rbind2(x0, object$x)
  return(x)
}

#' Predict method for `sglssnal` objects
#' @param object Fitted object of class `sglssnal`.
#' @param newdata A matrix of new data.
#' @param ... Additional arguments passed to or from other methods.
#' @return A numeric matrix of predictions, with columns corresponding to lambda values.
#' @method predict sglssnal
#' @export
predict.sglssnal <- function(object, newdata, ...) {
  result <- newdata %*% object$x
  if (object$intercept) {
    result <- result + object$x0
  }
  return(result)
}
