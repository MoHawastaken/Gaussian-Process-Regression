#'Simulation for Regression - Deterministic ground truth
#'
#'Simulates Regression problems with arbitrary regression function on a
#'hyperrectangle of arbitrary dimension and analyzes the quality of the Gaussian
#'process fit.
#'
#'When the regression function \code{func} is defined on a hyperrectangle in
#'\eqn{R^D}, limits has to be a \eqn{D x 2} matrix where the k-th line contains
#'the lower and upper bound of the hyperrectangle in the k-th dimension. If the
#'argument \code{training_points} is used, it should be a \eqn{D x n} matrix
#'with each column representing a training point. If \code{training_points} is
#'missing, \code{training_size} training points are drawn randomly from the
#'uniform distribution on the hyperrectangle determined by limits.
#'
#'The noise function given in \code{observation_noise} determines the noise
#'which is added to the values of func at the training points before Gaussian
#'process regression is done. When called with the matrix of training points,
#'\code{observation_noise} is supposed to return a vector of length n,
#'containing the error. Suitable iid. noise functions can be created from
#'standard distributions with \code{iid_noise}. See the examples for an example
#'of non iid. noise.
#'
#'The predictions of the Gaussian process are compared to the ground truth
#'(function values of func) on an equispaced grid of \code{test_size} points in
#'the hyperrectangle. Before returning the \code{summary} of the absolute error
#'in these test points, \code{simulate_regression} provides two plots. One of
#'the absolute error over the predicted variance and one of the regression
#'function \code{func} and the predictions of the Gaussian process. For some
#'inputs, the first plot may not be very informative. Setting \code{show_pred =
#'FALSE}, the first plot is not shown. If \eqn{D > 1}, the second plot is made
#'with respect to the first dimension, with all other variables fixed at the
#'middle of the their corresponding interval.
#'
#'@usage \preformatted{simulate_regression(func, limits, training_points,
#'  training_size = 10L, observation_noise = function(x) 0, test_size = 10000L,
#'  ...) }
#'@param func regression function
#'
#'@param limits matrix containing the limits of the hyperrectangle
#'
#'@param training_points matrix of points inside the hyperrectangle; optional
#'
#'@param training_size number of training points; unnecessary if training_points
#'  are given
#'
#'@param observation_noise function to generate the noise of the training data
#'
#'@param test_size number of test points at which regression function and
#'  Gaussian process are compared
#'
#'@param show_pred logical determining if plot of absolute prediction error over
#'  predicted variance is to be shown
#'
#'@param ... arguments passed on to the constructor of the Gaussian process
#'
#'@return \code{summary} of the absolute error of the Gaussian process
#'  predictions in the test points
#'
#'@examples
#'f <- function(x) 0.1 * x^3
#'limits <- matrix(c(-6, 6), nrow = 1)
#'X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
#'observation_noise <- iid_noise(rnorm, sd = 2)
#'simulate_regression(f, limits, X, observation_noise = observation_noise)
#'
#'f <- function(x) sum(x)^2
#'limits <- c(-1, 1, -1, 1, -1, 1)
#'observation_noise <- iid_noise(rcauchy)
#'simulate_regression(f, limits, training_size = 10, observation_noise = observation_noise)
#'
#'f <- function(x) 0.1 * sum(x^3)
#'limits <- matrix(c(-5, 5, -5, 5, -5, 5), ncol = 2, byrow = TRUE)
#'observation_noise <- function(X) {
#'  multivariate_normal(1, rep(0, ncol(X)), covariance_matrix(X, X, cov_func(sqrexp, l = 0.1)))
#'}
#'simulate_regression(f, limits, training_size = 10, observation_noise = observation_noise, noise = 1)
#'
#'@name simulate_regression
#'@export
simulate_regression <- function(func, limits, training_points, training_size = 10L, 
                    observation_noise = function(x) 0, test_size = 10000L, show_pred = TRUE, ...) {
  # Check correctness of inputs
  stopifnot(is.function(func), is.numeric(limits), length(limits) %% 2 == 0)
  stopifnot(is.numeric(training_size), training_size > 0)
  stopifnot(is.function(observation_noise), is.numeric(test_size), test_size > 0)
  if (!is.matrix(limits)) limits <- matrix(limits, ncol = 2, byrow = TRUE)
  
  # Training. Constructor uses fit to get best Gaussian model
  D <- nrow(limits)
  if (missing(training_points)) {
    training_points <- t(sapply(1:D, function(i) runif(training_size, limits[i,1], limits[i,2])))
  } else {
    stopifnot(nrow(limits) == nrow(training_points))
    stopifnot(all(training_points >= limits[, 1]), 
              all(training_points <= limits[, 2]))
  }
  y <- apply(training_points, 2, func) + observation_noise(training_points)
  Gaussian <- GPR$new(training_points, drop(y), ...)
  
  # Test the model on a large set of test points (size test_size)
  test_points <- combine_all(lapply(1:D, 
                  function(i) seq(limits[i, 1], limits[i, 2], length.out = test_size^(1/D))))
  predictions <- Gaussian$predict(test_points, pointwise_var = TRUE)
  residual <- predictions[, 1] - apply(test_points, 2, func)
  cat("The mean absolute difference of predictions and ground truth",
              "in the considered limits is ", mean(abs(residual)), "\n")
  
  # Plot error over predicted variance
  if (show_pred){
    variance <- predictions[, 2]
    plot(variance, abs(residual), xlim = c(0, 1.1 * variance[length(variance)]),
         main = "Connection of absolute prediction error \n and predicted variance",
         xlab = "Variance", ylab = "Absolute Prediction Error", col = "blue")
  }
  # Plot regression function and estimated function.
  x <- seq(limits[1, 1], limits[1, 2], by = 0.05)
  if (D == 1) {
    predictions <- Gaussian$predict(x, pointwise_var = TRUE)
    ground_truth <- sapply(x, func)
  } else {
    message("Plot with respect to first variable with all others fixed.")
    plot_points <- rbind(x, matrix(apply(limits, 1, mean)[-1], 
                                          nrow = D - 1, ncol = length(x)))
    ground_truth <- apply(plot_points, 2, func)
    predictions <- Gaussian$predict(plot_points, pointwise_var = TRUE)
  }
  dat <- data.frame(x = x, ground_truth = ground_truth, 
                    regression = predictions[, 1], variance = predictions[, 2])
  dat1 <- tidyr::gather(dat, ground_truth, regression, key = "Function", value = "value")
  print(ggplot2::ggplot(dat1, ggplot2::aes(x = x, y = value, colour = Function)) +
    ggplot2::theme_classic() +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(inherit.aes = F, 
        data = dat,
        mapping = ggplot2::aes(x = x, ymin = regression - 2*sqrt(pmax(variance,0)),
                              ymax = regression + 2*sqrt(pmax(variance, 0))), 
        alpha = 0.3))
  return(summary(abs(residual)))
}

#'Simulation for Regression - Random ground truth
#'
#'Simulates Regression problems on a hyperrectangle of arbitrary dimension with
#'the ground truth drawn from a Gaussian process and analyzes the quality of the
#'Gaussian process fit.
#'
#'\code{limits} should be a \eqn{D x 2} matrix to model a Gaussian process
#'in \eqn{\mathbb{R}^D}, the k-th line containing the lower and upper bound of
#'the hyperrectangle in the k-th dimension. The ground truth value at the test
#'points (an equispaced grid of \code{test_size} points in the hyperrectangle)
#'is randomly drawn from a multivariate distribution with mean 0 and the
#'covariance matrix determined by the covariance function \code{actual_cov}.
#'Then \code{training_size} training points are drawn randomly from this grid of
#'test points.
#'
#'The noise function given in \code{observation_noise} determines the noise
#'which is added to the ground truth values at the training points before
#'Gaussian process regression is done. When called with the matrix of training
#'points, \code{observation_noise} is supposed to return a vector of length n,
#'containing the error. Suitable iid. noise functions can be created from
#'standard distributions with \code{iid_noise}. See the examples for an example
#'of non iid. noise.
#'
#'The predictions of the Gaussian process at the test points are compared to the
#'ground truth. Before returning the \code{summary} of the absolute error in
#'these test points, \code{simulate_regression_gp} provides up to two plots. If
#'\code{show_pred = TRUE}, one plot shows the absolute error over the predicted
#'variance and, if \eqn{D = 1} another shows the ground truth GP and the
#'predictions of the fitted GP.
#'
#'@usage \preformatted{simulate_regression_gp(actual_cov, limits,
#'  observation_noise = function(x) 0, test_size = 300, training_size = 10,
#'  random_training = TRUE, show_pred = FALSE, ...) }
#'@param actual_cov covariance function with which ground truth GP is generated
#'
#'@param limits matrix containing the limits of the hyperrectangle
#'
#'@param observation_noise function to generate the noise of the training data
#'
#'@param test_size integer, number of test points at which ground truth and
#'  fitted Gaussian process are compared
#'
#'@param training_size integer, number of training points
#'
#'@param random_training logical; if TRUEtraining points are to be drawn from
#'  the uniform distribution on the hyperrectangle instead of being chosen
#'  equispaced
#'
#'@param show_pred logical; if TRUE the plot of the absolute prediction error
#'  over the predicted variance is to be shown
#'
#'@param ... arguments passed on to the constructor of the Gaussian process
#'
#'@return \code{summary} of the absolute error of the Gaussian process
#'  predictions in the test points
#'
#'@examples
#'simulate_regression_gp(cov_func(sqrexp, l = 1), limits = c(-5, 5), training_size = 10)
#'
#'simulate_regression_gp(cov_func(rationalquadratic, alpha = 0.5, l = 1), 
#'  limits = c(-10, 10), observation_noise = iid_noise(rnorm, 0.2), 
#'   training_size = 50, noise = 0.2)
#'
#'observation_noise <- function(X) {
#'  multivariate_normal(1, rep(0, ncol(X)),
#'         covariance_matrix(X, X, function(x,y) sqrexp(x, y, 0.1)))
#'}
#'simulate_regression_gp(cov_func(sqrexp, l = 1), limits = c(-5, 5), training_size = 10,
#'  observation_noise = observation_noise)
#'@name simulate_regression_gp
#'@export
simulate_regression_gp <- function(actual_cov, limits, observation_noise = function(x) 0, test_size = 300L, 
                training_size = 10L, random_training = TRUE, show_pred = FALSE, ...) {
  # Check correctness of inputs
  stopifnot(is.function(actual_cov), is.numeric(limits), length(limits) %% 2 == 0)
  stopifnot(is.function(observation_noise), is.numeric(test_size), test_size > 0)
  stopifnot(is.numeric(training_size), training_size > 0, training_size < test_size)
  stopifnot(is.logical(random_training))
  if (!is.matrix(limits)) limits <- matrix(limits, ncol = 2, byrow = TRUE)
  
  # Draw Gaussian process f with covariance function actual_cov
  D <- nrow(limits)
  testpoints <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], 
                                                         length.out = test_size^(1/D))))
  K <- covariance_matrix(testpoints, testpoints, actual_cov)
  f <- multivariate_normal(1, rep(0, nrow(K)), K)
  
  # Training
  if (random_training) training_set <- sample.int(ncol(testpoints), training_size)
  else training_set <- (1:training_size)*floor(ncol(testpoints)/training_size)
  X <- testpoints[, training_set, drop = FALSE]
  y <- f[training_set] + observation_noise(X)
  regression_GP <- GPR$new(X, drop(y), ...)
  
  # Testing
  if (D == 1) {
    plot_l <- regression_GP$plot(drop(limits), test_size)
    residual <- plot_l$pred[, 1] - f
    variance <- plot_l$pred[, 2]
  } else {
    prediction <- regression_GP$predict(testpoints)
    residual <- prediction[, 1] - f
    variance <- prediction[, 2]
  }
  if (show_pred){
    plot(variance, abs(residual), xlim = c(0, 1.1 * variance[length(variance)]),
         main = "Connection of absolute prediction error \n and predicted variance",
         xlab = "Variance", ylab = "Absolute Prediction Error", col = "blue")
  }
  if (D == 1) print(plot_l$plot + ggplot2::geom_line(data = data.frame(x = drop(testpoints), y = f),
                                                     ggplot2::aes(x = x, y = y), colour = "green"))
  message(paste("The mean absolute difference of predictions and ground truth",
      "in the considered limits is ", mean(abs(residual)), "\n"))
  return(summary(abs(residual)))
}

#'Simulation for Classification
#'
#'Simulates Classification problems with arbitrary classification function on a
#'hyperrectangle of arbitrary dimension and analyzes the quality of the Gaussian
#'process predictions.
#'
#'When the classification function \code{func} is defined on a hyperrectangle in
#'\eqn{R^D}, limits has to be a \eqn{D x 2} matrix where the k-th line contains
#'the lower and upper bound of the hyperrectangle in the k-th dimension. If the
#'argument \code{training_points} is used, it should be a \eqn{D x n} matrix
#'with each column representing a training point. If \code{training_points} is
#'missing, \code{training_size} training points are drawn randomly from the
#'uniform distribution on the hyperrectangle determined by limits.
#'
#'The predictions of the Gaussian process are compared to the ground truth
#'(function values of func) on an equispaced grid of \code{test_size} points in
#'the hyperrectangle. If \eqn{D <= 2}, a plot of the decision regions of the
#'Gaussian Process along with the true classification of the testpoints is
#'provided. If \eqn{D = 1}, the predicted posterior probability is also plotted.
#'
#'@usage \preformatted{simulate_classification <- function(func, limits,
#'  training_points, training_size = 10L, test_size = 10000L, ...)}
#'
#'@param func A function with values -1 and 1 defined on the hyperrectangle given
#'  by limits.
#'
#'@param limits matrix containing the limits of the hyperrectangle
#'
#'@param training_points matrix of points inside the hyperrectangle; optional
#'
#'@param training_size number of training points; unnecessary if training_points
#'  are given
#'
#'@param test_size number of test points at which classification function and
#'  Gaussian process predictions are compared
#'
#'@param ... arguments passed on to the constructor of the Gaussian process
#'
#'@return \code{summary} of the absolute error of the Gaussian process
#'  predictions in the test points
#'
#' @examples
#'  f <- function(x) (x < - 2) + (x > 2) - (-2 <= x && x <= 2)
#'  limits <- matrix(c(-4, 4), nrow = 1, byrow = TRUE)
#'  k <- function(x, y) sqrexp(x, y, 1)
#'  simulate_classification(func = f, limits = limits, training_size = 20, k = k)
#'
#'  f <- function(x) (sum(abs(x)) > 2.5) - (!(sum(abs(x)) > 2.5))
#'  limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
#'  k <- function(x, y) sqrexp(x, y, 1)
#'  simulate_classification(func = f, limits = limits, training_size = 50, k = k)
#'  
#'@name simulate_classification
#'@export
simulate_classification <- function(func, limits, training_points,
                                    training_size = 10L, test_size = 10000L, ...) {
  # Check correctness of inputs
  stopifnot(is.function(func), is.numeric(limits), length(limits) %% 2 == 0)
  stopifnot(is.numeric(training_size), training_size > 0)
  if (!is.matrix(limits)) limits <- matrix(limits, ncol = 2, byrow = TRUE)
  
  # Training
  D <- nrow(limits)
  if (missing(training_points)) {
    training_points <- t(sapply(1:nrow(limits), 
                        function(i) runif(training_size, limits[i, 1], limits[i, 2])))
  } else stopifnot(nrow(limits) == nrow(training_points))
  y <- apply(training_points, 2, func)
  Gaussian <- GPC$new(training_points, y, ...)
  
  # Testing
  test_points <- combine_all(lapply(1:D, 
                        function(i) seq(limits[i, 1], limits[i, 2], length.out = test_size^(1/D))))
  plot_l <- Gaussian$plot(testpoints = test_points)
  predictions <- plot_l$pred
  residual <- 2 * as.integer(predictions >= 0.5) - 1 - apply(test_points, 2, func)
  message(paste("Proportion of misclassified test points.", mean(abs(residual))/2, "\n"))
  print(plot_l$plot)
  return(summary(abs(residual)))
}
 
combine_all <- function(lst) {
  l <- length(lst)
  lengths <- sapply(lst, length)
  rev_lengths <- lengths[l:1]
  prods <- c(1, cumprod(lengths))
  rev_prods <- c(1, cumprod(rev_lengths))
  out <- matrix(0, nrow = l, ncol = prods[l+1])
  for (k in 1:l) {
    out[k, ] <- rep(lst[[k]], each = rev_prods[l + 1 - k], times = prods[k])
  }
  out
}

#'Independent and identically distributed noise
#'
#'Creates suitable noise functions for \code{simulate_regression} and
#'\code{simulate_regression_gp} with iid noise generated by standard
#'distributions.
#'
#'The first argument of \code{distribution} has to be the number of random
#'numbers to generate.
#'
#'@usage \preformatted{iid_noise(distribution, ...)}
#'
#'@param distribution A function generating random numbers

#'@param ... arguments passed on to \code{distribution}
#'
#'@return a suitable noise function for \code{simulate_regression} and
#'\code{simulate_regression_gp}
#'
#' @examples
#'  observation_noise <- iid(rnorm, sd = 0.1)
#'  observation_noise <- iid(rcauchy)
#'
#'@name iid_noise
#'@export
iid_noise <- function(distribution, ...) {
  force(distribution)
  function(X) distribution(ncol(X), ...)
}

#'Example plots of Gaussian Process simulations
#'
#'Runs examples of \code{simulate_regression}, \code{simulate_regression_gp} and
#'\code{simulate_classification}
#'
#'@usage \preformatted{examples()}

#'@return Nothing
#'@name examples
#' @export
examples <- function() {
  # Regression
  # example 1: one-dimensional
  f <- function(x) 0.1 * x^3
  limits <- matrix(c(-6, 6), nrow = 1)
  X <- matrix(seq(-5, 5, by = 0.2), nrow = 1)
  observation_noise <- iid_noise(rnorm, sd = 2)
  simulate_regression(f, limits, X, noise = 1, observation_noise = observation_noise)
  
  # example 2: cauchy distributed noise
  f <- function(x) sin(10 * x)
  limits <- matrix(c(0, 1), nrow = 1)
  X <- matrix(seq(0, 1, by = 0.05), nrow = 1)
  observation_noise <- iid_noise(rcauchy)
  simulate_regression(f, limits, X, noise = 1, observation_noise = observation_noise)
  
  # example 3: two-dimensional
  f <- function(x) 0.1 * sum(x^2)
  limits <- matrix(c(-5.5, 5.5, -5.5, 5.5), nrow = 2, byrow = TRUE)
  X <- combine_all(list(seq(-5, 5, by = 1), seq(-5 ,5, by = 1)))
  observation_noise <- iid_noise(rnorm, sd = 1)
  simulate_regression(f, limits, X, noise = 1, observation_noise = observation_noise)
  
  # Regression for Gaussian processes
  simulate_regression_gp(cov_func(sqrexp, l = 1), limits = matrix(c(-5, 5), nrow = 1), 
                         training_size = 10, random_training = TRUE, noise = 0.1)
  
  simulate_regression_gp(cov_func(sqrexp, l = 0.1), limits = matrix(c(-5, 5), nrow = 1), 
                         training_size = 10, random_training = TRUE, noise = 0.1)
  
  simulate_regression_gp(cov_func(sqrexp, l = 1), limits = matrix(c(-5, 5), nrow = 1), 
                         training_size = 10, random_training = TRUE, noise = 1)
  
  simulate_regression_gp(cov_func(polynomial, sigma = 1, p = 3), limits = matrix(c(-5, 5), nrow = 1), 
                         training_size = 10, random_training = TRUE, noise = 1)
  
  # Classification
  # example 1
  f <- function(x) (x < - 2) + (x > 2) - (-2 <= x && x <= 2)
  limits <- matrix(c(-4, 4), nrow = 1, byrow = TRUE)
  k <- cov_func(sqrexp, l = 1)
  simulate_classification(func = f, limits = limits, training_size = 20, k = k)
  
  # example 2
  f <- function(x) (sum(abs(x)) > 2.5) - (!(sum(abs(x)) > 2.5))
  limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
  k <-  cov_func(sqrexp, l = 1)
  simulate_classification(func = f, limits = limits, training_size = 50, k = k)
  
  # Haus
  f <- function(x) {
    ((sum(abs(x)) < 4 && x[2] > 1)|| (x[1] > -3 && x[1] < 3 && x[2] <= 1)) - 
      (!((sum(abs(x)) < 4 && x[2] > 1)|| (x[1] > -3 && x[1] < 3 && x[2] <= 1))) -
      2*(x[1] > -0.75 && x[1] < 0.75 && x[2] < -2) -
      2*(x[1] > -2 && x[1] < -0.5 && x[2] > -1 && x[2] < 0.5) -
      2*(x[1] > 0.5 && x[1] < 2 && x[2] > -1 && x[2] < 0.5)
  }
  limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
  k <-  cov_func(sqrexp, l = 0.5)
  simulate_classification(func = f, limits = limits, training_size = 600, k = k)
  
  # Example R
  g <- function(x){
    if (x == 0) return(-1)
    if (x == -3) return(-1)
    return(x)
  }
  f <- function(x) {
    g(((x[1] > -3.2 && x[1] < 1.8 && x[2] <= 3.2 && x[2] > 0) || (x[1] > -3 && x[1] < 3.5 && x[2] <= 3 && x[2] <= 0) ) - 
        (!((x[1] > -3.2 && x[1] < 1.8 && x[2] <= 3.2 && x[2] > 0) || (x[1] > -3 && x[1] < 3.5 && x[2] <= 3 && x[2] <= 0) )) -
        2*(x[1] > -1.5 && x[1] < 0 && x[2] < -1) -
        2*(x[1] > 0 && x[1] < 2 && x[2]  < -1 && sum(x) < -1.8) -
        2*(x[1] > -2 && x[1] < 0.5 && x[2] > 0 && x[2] < 2) -
        2*(x[1] > 1 && x[2] < 0 && x[2] > -3.5 && x[1] < 3.5 && sum(x) > -0.4))
  }
  
  limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
  k <-  cov_func(sqrexp, l = 1)
  simulate_classification(func = f, limits = limits, training_size = 600, k = k)
}
