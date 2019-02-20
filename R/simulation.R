#limits D x 2 Matrix, X D x N
max_predict <- 10000

#' Simulation for Regression
#' 
#' Simulates data for Regression problems, which can then be analyzed via a Gaussian process
#' 
#' @usage \preformatted{simulate_regression(func, limits, training_points, num_data = 10, 
#' noise = 0, error = function(x) 0, cov_names = names(cov_dict))
#'}
#' @param func
#' 
#' @examples 
#' 
#' @name simulate_regression
#' @export
simulate_regression <- function(func, limits, training_points, num_data = 10, 
                                noise = 0, error = function(x) 0, cov_names = names(cov_dict)) {
  D <- nrow(limits)
  if (missing(training_points)) {
    training_points <- t(sapply(1:D, function(i) runif(num_data, limits[i,1], limits[i,2])))
  } else stopifnot(nrow(limits) == nrow(training_points))
  
  # Use fit to get best Gaussian model
  y <- apply(training_points, 2, func) + error(ncol(training_points))
  Gaussian <- GPR$new(training_points, y, noise = noise, cov_names = cov_names)
  
  # Test the model on a large set of test points (size max_predict)
  test_points <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], length.out = max_predict^(1/D))))
  predictions <- Gaussian$predict(test_points, pointwise_var = TRUE)
  residual <- predictions[, 1] - apply(test_points, 2, func)
  cat("The mean absolute difference of predictions and ground truth",
              "in the considered limits is ", mean(abs(residual)), "\n")
  
  # Plot error over predicted variance
  variance <- predictions[, 2]
  resid <- abs(residual)
  plot(variance, resid, xlim = c(0, 1.1*variance[length(variance)]),
       main = "Connection of absolute prediction error \n and predicted variance",
       xlab = "Variance", ylab = "Absolute Prediction Error", col = "blue")
  
  # Plot regression function and estimated function.
  x <- seq(limits[1, 1], limits[1, 2], by = 0.05)
  if (D == 1) {
    predictions <- Gaussian$predict(x, pointwise_var = TRUE)
    ground_truth <- sapply(x, func)
  } else {
    print("Plot with respect to first variable with all others fixed.")
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
}

max_set <- 300
#' @export
simulate_regression_gp <- function(actual_cov, noise, limits, n_training_data = 10, 
                                   random_training = TRUE, regression_noise = noise,...) {
  if (!is.matrix(limits)) dim(limits) <- c(length(limits) / 2, 2)
  D <- nrow(limits)
  testpoints <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], 
                                                         length.out = max_set^(1/D))))
  K <- covariance_matrix(testpoints, testpoints, actual_cov)
  f <- multivariate_normal(1, rep(0, nrow(K)), K)
  actual_GP <- GPR$new(testpoints, drop(f), noise = 0.1, k = actual_cov)
  
  if (random_training) X <- t(sapply(1:D, function(i) runif(n_training_data, limits[i,1], limits[i,2])))
  else X <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], 
                                                    length.out = n_training_data^(1/D))))
  y <- actual_GP$predict(X)[, 1] + rnorm(n_training_data, 0, noise)
  regression_GP <- GPR$new(X, y, noise = regression_noise, ...)
  testpoints <- drop(testpoints)
  lst <- regression_GP$plot(testpoints)
  lst$plot + ggplot2::geom_line(data = data.frame(x = testpoints, y = f),
                         ggplot2::aes(x = x, y = y), colour = "green")
}

#' Simulation for Classification
#' 
#' Simulates data for Classification problems, which can then be analyzed via a Gaussian process
#' 
#' @usage \preformatted{simulate_classification(func, training_points, limits, k, num_data = 10)
#'}
#' @param func A function for the simulation of data points
#' @param training_points A vector of training points
#' @param limits 
#' 
#' @examples 
#' 
#' @name simulate_classification
#' @export
simulate_classification <- function(func, training_points, limits, k, num_data = 10) {
  D <- nrow(limits)
  if (missing(training_points)) {
    training_points <- t(sapply(1:nrow(limits), function(i) runif(num_data, limits[i,1], limits[i,2])))
  } else stopifnot(nrow(limits) == nrow(training_points))
  test_points <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], length.out = max_predict^(1/D))))
  y <- apply(training_points, 2, func)
  Gaussian <- GPC$new(training_points, y, 1e-5, k)
  plot_l <- Gaussian$plot(test_points)
  predictions <- plot_l$pred
  residual <- 2*as.integer(predictions >= 0.5) - 1 - apply(test_points, 2, func)
  message(paste("Proportion of misclassified test points.", mean(abs(residual))/2, "\n"))
  print(plot_l$plot)
}
 
#' @export
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

#' @export
normal <- function(sd, mean = 0) {
  function(k) rnorm(k, mean = mean, sd = sd)
}

#' @export
exp_distr <- function(rate) {
  function(k) rexp(k, rate)
}

#' @export
beta_distr <- function(shape1, shape2, ncp = 0){
  function(k) rbeta(k, shape1, shape2, ncp)
}

#' @export
examples <- function() {
  # Regression
  # example 1: one-dimensional
  f <- function(x) 0.1*x^3
  limits <- matrix(c(-6, 6), nrow = 1)
  X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
  error <- normal(2)
  simulate_regression(f, limits, X, noise = 1, error = error)
  
  # example 2: beta distributed error
  f <- function(x) sin(10*x)
  limits <- matrix(c(0, 1), nrow = 1)
  X <- matrix(seq(0,1,by = 0.05), nrow = 1)
  error <- beta_distr(2,3)
  simulate_regression(f, limits, X, noise = 1, error = error)
  
  # example 3: two-dimensional
  f <- function(x) 0.1*sum(x^2)
  limits <- matrix(c(-5.5, 5.5, -5.5, 5.5), nrow = 2, byrow = TRUE)
  X <- combine_all(list(seq(-5,5,by = 1), seq(-5,5,by = 1)))
  error <- normal(1)
  simulate_regression(f, limits, X, noise = 1, error = error)

  
  # Classification
  # example 1
  f <- function(x) (x < - 2) + (x > 2) - (-2 <= x && x <= 2)
  limits <- matrix(c(-4, 4), nrow = 1, byrow = TRUE)
  k <- function(x, y) sqrexp(x, y, 1)
  simulate_classification(func = f, limits = limits, k = k, num_data = 20)
  
  # example 2
  f <- function(x) (sum(abs(x)) > 2.5) - (!(sum(abs(x)) > 2.5))
  limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
  k <- function(x, y) sqrexp(x, y, 1)
  simulate_classification(func = f, limits = limits, k = k, num_data = 50)
  
  # Haus
  f <- function(x) {
    ((sum(abs(x)) < 4 && x[2] > 1)|| (x[1] > -3 && x[1] < 3 && x[2] <= 1)) - 
      (!((sum(abs(x)) < 4 && x[2] > 1)|| (x[1] > -3 && x[1] < 3 && x[2] <= 1))) -
      2*(x[1] > -0.75 && x[1] < 0.75 && x[2] < -2) -
      2*(x[1] > -2 && x[1] < -0.5 && x[2] > -1 && x[2] < 0.5) -
      2*(x[1] > 0.5 && x[1] < 2 && x[2] > -1 && x[2] < 0.5)
  }
  limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
  k <- function(x, y) sqrexp(x, y, 1)
  simulate_classification(func = f, limits = limits, k = k, num_data = 600)
  
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
  k <- function(x, y) sqrexp(x, y, 1)
  simulate_classification(func = f, limits = limits, k = k, num_data = 600)
}
