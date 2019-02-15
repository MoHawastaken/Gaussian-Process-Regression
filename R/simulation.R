 multivariate_normal <- function(n, mean, covariance) {
  stopifnot(is.numeric(mean), is.numeric(covariance), length(mean) == nrow(covariance))
  L <- t(chol(covariance))
  c(mean) + L%*%matrix(rnorm(n*length(mean), 0, 1), nrow = length(mean))
 }
 
 #limits D x 2 Matrix, X D x N
 simulate <- function(func, limits, training_points, noise = 0, error = function(x) 0, cov_names = names(cov_dict)) {
   y <- apply(training_points, 2, func) + error(ncol(training_points))
   best <- fit(training_points, y, noise, cov_names)
   k <- function(x,y) do.call(usedcov$func[[1]], append(list(x, y), v))
   #Gaussian <- GPR$
 }
 
normal <- function(sd, mean = 0) {
  function(k) rnorm(k, mean = mean, sd = sd)
}
error <- normal(0.1)