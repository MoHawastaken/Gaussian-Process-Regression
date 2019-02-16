multivariate_normal <- function(n, mean, covariance) {
  stopifnot(is.numeric(mean), is.numeric(covariance), length(mean) == nrow(covariance))
  L <- t(chol(covariance))
  c(mean) + L%*%matrix(rnorm(n*length(mean), 0, 1), nrow = length(mean))
}
 
#limits D x 2 Matrix, X D x N
max_predict <- 10000
simulate <- function(func, training_points, limits, noise = 0, error = function(x) 0, cov_names = names(cov_dict)) {
  stopifnot(nrow(limits) == nrow(training_points))
  D <- nrow(limits)
  test_points <- combine_all(lapply(1:D, function(i) seq(limits[i, 1], limits[i, 2], length.out = max_predict^(1/D))))
  y <- apply(training_points, 2, func) + error(ncol(training_points))
  print(1)
  best <- fit(training_points, y, noise, cov_names)
  print(2)
  k <- function(x,y) do.call(cov_df[best$cov, ]$func[[1]], append(list(x, y), best$par))
  Gaussian <- GPR$new(X, y, k, noise)
  predictions <- Gaussian$predict(test_points, pointwise_var = TRUE)
  residual <- predictions[, 1] - apply(test_points, 2, func)
  cat("The mean absolute difference of predictions and ground truth",
              "in the considered limits is ", mean(abs(residual)), "\n")
  plot(predictions[, 2], abs(residual))
  x <- seq(limits[1, 1], limits[1, 2], by = 0.05)
  if (D == 1) y <- cbind(sapply(x, func), Gaussian$predict(x))
  else {
    print("Plot with respect to first variable with all others fixed.")
    test_points <- rbind(test_points, matrix(apply(limits, 1, mean)[-1], nrow = D-1))
    y <- cbind(apply(test_points, 2, func), Gaussian$predict(test_points, pointwise_var = TRUE))
  }
  print(y)
  dat <- data.frame(x = x, y = y)
  dat1 <- tidyr::gather(dat, y.1, y.2, key = "Function", value = "value")
  ggplot2::ggplot(dat1, ggplot2::aes(x = x, y = value, colour = Function)) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous("output, f(x)") +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(inherit.aes = F, mapping = ggplot2::aes(x = x, ymin = y.2 - 2*sqrt(pmax(y.3,0)),
                                ymax = y.2 + 2*sqrt(pmax(y.3,0))), data = dat, alpha = 0.3) 
  
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
 
normal <- function(sd, mean = 0) {
  function(k) rnorm(k, mean = mean, sd = sd)
}

# example 1
f <- function(x) 0.1*x^3
limits <- matrix(c(-5.5, 5.5), nrow = 1)
X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
error <- normal(1)
simulate(f, X, limits, noise = 1, error = error, cov_names = c("gammaexp", "rationalquadratic"))

# example 1
f <- function(x) 0.1*sum(x^2)
limits <- matrix(c(-5.5, 5.5, -5.5, 5.5), nrow = 2, byrow = TRUE)
X <- combine_all(list(seq(-5,5,by = 0.2), seq(-5,5,by = 0.2)))
error <- normal(1)
#Bei Fit ergibt sich ein error
#simulate(f, X, limits, noise = 1, error = error, cov_names = c("gammaexp", "rationalquadratic"))
