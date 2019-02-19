context("test-gpc")

test_that("gpc works", {
  # most basic classification
  X <- matrix(seq(-1,1,by = 0.1), nrow = 1)
  y <- 2*as.integer(X > 0) - 1
  kappa <- function(x,y) exp(-3*(x - y)^2)
  gaussian_classifier <- GPC$new(X, y, 1e-5, kappa)
  expect_lt(gaussian_classifier$predict_class(-0.2), 0.5)
  expect_gt(gaussian_classifier$predict_class(0.2), 0.5)
  
  # more density for negatives
  X <- matrix(c(seq(-1, -0.1, by = 0.1), seq(0, 1, by = 0.2)), nrow = 1)
  y <- 2*as.integer(X > 0) - 1
  kappa <- function(x,y) exp(-3*(x - y)^2)
  gaussian_classifier <- GPC$new(X, y, 1e-5, kappa)
  expect_lt(gaussian_classifier$predict_class(-0.2), 0.5)
  expect_gt(gaussian_classifier$predict_class(0.2), 0.5)
  
  #2-dim Raster
  s <- seq(-1, 1, by = 0.5)
  X <- matrix(c(rep(s, each = length(s)), rep(s, times = length(s))), nrow = 2, byrow = T)
  y <- 2*as.integer(X[1, ] > X[2, ]) - 1
  kappa <- function(x,y) sqrexp(x,y,l=1)
  gaussian_classifier <- GPC$new(X, y, 1e-5, kappa)
  expect_lt(gaussian_classifier$predict_class(matrix(c(0,1), nrow = 2)), 0.5)
  expect_gt(gaussian_classifier$predict_class(matrix(c(-0.3,-0.9), nrow = 2)), 0.5)
  
  #2-dim normal distr. clusters
  n <- 10
  X <- cbind(multivariate_normal(n, c(0.5,0.5), diag(c(0.1,0.1))), multivariate_normal(n, c(-0.5,-0.5), diag(c(0.1,0.1))))
  y <- rep(c(1,-1), each = n)
  kappa <- function(x,y) sqrexp(x,y,l=1)
  gaussian_classifier <- GPC$new(X, y, 1e-5, kappa)
  expect_lt(gaussian_classifier$predict_class(matrix(c(-0.2,-0.2), nrow = 2)), 0.5)
  expect_gt(gaussian_classifier$predict_class(matrix(c(0.2,0.2), nrow = 2)), 0.5)
})
