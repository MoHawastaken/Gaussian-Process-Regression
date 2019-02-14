multivariate_normal <- function(n, mean, covariance) {
  stopifnot(is.numeric(mean), is.numeric(covariance), length(mean) == nrow(covariance))
  L <- t(chol(covariance))
  c(mean) + L%*%matrix(rnorm(n*length(mean), 0, 1), nrow = length(mean))
}