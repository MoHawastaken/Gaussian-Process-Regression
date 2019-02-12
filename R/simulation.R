multivariate_normal <- function(n, mean, covariance) {
  stopifnot(is.numeric(mean), is.numeric(covariance), length(mean) == nrow(covariance))
  L <- t(chol(covariance))
  mean + L%*%rnorm(n, 0, 1)
}