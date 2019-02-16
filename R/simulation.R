#limits D x 2 Matrix, X D x N
max_predict <- 10000
#' @export
simulate_regression <- function(func, limits, training_points, num_data = 10, 
                                noise = 0, error = function(x) 0, cov_names = names(cov_dict)) {
  D <- nrow(limits)
  if (missing(training_points)) {
    training_points <- t(sapply(1:D, function(i) runif(num_data, limits[i,1], limits[i,2])))
  } else stopifnot(nrow(limits) == nrow(training_points))
  test_points <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], length.out = max_predict^(1/D))))
  y <- apply(training_points, 2, func) + error(ncol(training_points))
  print(1)
  best <- fit(training_points, y, noise, cov_names)
  print(2)
  k <- function(x,y) do.call(cov_df[best$cov, ]$func[[1]], append(list(x, y), best$par))
  Gaussian <- GPR$new(training_points, y, k, noise)
  predictions <- Gaussian$predict(test_points, pointwise_var = TRUE)
  residual <- predictions[, 1] - apply(test_points, 2, func)
  # Visualizations
  cat("The mean absolute difference of predictions and ground truth",
              "in the considered limits is ", mean(abs(residual)), "\n")
  plot(predictions[, 2], abs(residual))
  x <- seq(limits[1, 1], limits[1, 2], by = 0.05)
  if (D == 1) {
    predictions <- Gaussian$predict(x, pointwise_var = TRUE)
    ground_truth <- sapply(x, func)
  } else {
    print("Plot with respect to first variable with all others fixed.")
    plot_points <- rbind(x, matrix(apply(limits, 1, mean)[-1], 
                                          nrow = D-1, ncol = length(x)))
    ground_truth <- apply(plot_points, 2, func)
    predictions <- Gaussian$predict(plot_points, pointwise_var = TRUE)
  }
  dat <- data.frame(x = x, ground_truth = ground_truth, 
                    regression = predictions[, 1], variance = predictions[, 2])
  dat1 <- tidyr::gather(dat, ground_truth, regression, key = "Function", value = "value")
  ggplot2::ggplot(dat1, ggplot2::aes(x = x, y = value, colour = Function)) +
    ggplot2::theme_classic() +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(inherit.aes = F, 
        data = dat,
        mapping = ggplot2::aes(x = x, ymin = regression - 2*sqrt(pmax(variance,0)),
                              ymax = regression + 2*sqrt(pmax(variance,0))), 
        alpha = 0.3) 
}

<<<<<<< HEAD
#simulate_classification 
=======


simulate_classification <- function(func, training_points, limits, k, num_data = 10) {
  D <- nrow(limits)
  if (missing(training_points)) {
    training_points <- t(sapply(1:nrow(limits), function(i) runif(num_data, limits[i,1], limits[i,2])))
  } else stopifnot(nrow(limits) == nrow(training_points))
  test_points <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], length.out = max_predict^(1/D))))
  y <- apply(training_points, 2, func)
  Gaussian <- GPC$new(training_points, y, k, 1e-5)
  predictions <- Gaussian$predict_class2(test_points)
  residual <- 2*as.integer(predictions >= 0.5) - 1 - apply(test_points, 2, func)
  # Visualizations
  cat("Proportion of misclassified test points.", mean(abs(residual))/2, "\n")
  Gaussian$plot(test_points)
}
>>>>>>> 4b06cc0b86d0c64cd258605b65df599fe2e20eca
 
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

#Regression
# example 1
f <- function(x) 0.1*x^3
limits <- matrix(c(-6, 6), nrow = 1)
X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
error <- normal(1)
simulate_regression(f, X, limits, noise = 1, error = error, cov_names = c("gammaexp", "rationalquadratic"))

# example 1
f <- function(x) 0.1*sum(x^2)
limits <- matrix(c(-5.5, 5.5, -5.5, 5.5), nrow = 2, byrow = TRUE)
X <- combine_all(list(seq(-5,5,by = 1), seq(-5,5,by = 1)))
error <- normal(1)
simulate_regression(f, X, limits, noise = 1, error = error, cov_names = c("gammaexp", "rationalquadratic"))

#Classification
#example 1
f <- function(x) (x < - 2) + (x > 2) - (-2 <= x && x <= 2)
limits <- matrix(c(-4, 4), nrow = 1, byrow = TRUE)
k <- function(x, y) sqrexp(x, y, 1)
simulate_classification(func = f, limits = limits, k = k, num_data = 20)

#example 2
f <- function(x) (sum(abs(x)) > 2.5) - (!(sum(abs(x)) > 2.5))
limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
k <- function(x, y) sqrexp(x, y, 1)
simulate_classification(func = f, limits = limits, k = k, num_data = 50)

#Haus
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
