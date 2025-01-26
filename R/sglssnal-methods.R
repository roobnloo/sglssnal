#' Extract estimated coefficients from a `sglssnal` object
#' @param object Fitted object of class `sglssnal`.
#' @return A numeric vector of estimated coefficients.
#' @method coef sglssnal
#' @export
coef.sglssnal <- function(object) {
  object$x
}

#' Predict method for `sglssnal` objects
#' @param object Fitted object of class `sglssnal`.
#' @param newdata A matrix of new data.
#' @return A numeric vector of predictions.
#' @method predict sglssnal
#' @export
predict.sglssnal <- function(object, newdata) {
  newdata %*% stats::coef(object)
}
