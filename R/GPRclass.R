#'  Predictions and Plots for Gaussian process regression
#'
#'  Implements a gaussian process and gives tools to predict and plot its values for given testpoints
#' 
#'
#' @section Usage: 
#' \preformatted{GPR <- GPR$new(X, y, cov_Fun, noise)
#'
#'
#' GPR$predict(X*)
#' GPR$plot(testpoints)
#'}
#' @section Arguments:
#' 
#'   \code{X} matrix of inputs
#'
#'   \code{y} numeric vector of targets
#' 
#'   \code{cov_Fun} the predicted covarianz function of the gaussian process
#' 
#'   \code{noise} the predicted noise of the observations
#' 
#'   \code{X*} a numeric vector as the test input
#' 
#'   \code{testpoints} a matrix of testpoints
#'   
#'
#' @section Methods:
#' \code{$predict()} returns a numeric vector of the expected value of the underlying function f and their variance for the test input
#' 
#' \code{$plot()} displays the results of the predict function for all testpoints in a nice plot
#' 
#' 
#' @section Methods:
#' GPR has several subclasses where a covarianz function k(x,y) is given. The following subclasses are implemented:
#' 
#' \code{GPR <- GPR.constant$new(X, y, c, noise)} with \code{k(x,y) = c}
#' 
#' \code{GPR <- GPR.linear$new(X, y, cov_Fun, noise)} with \code{k(x,y) = sum(sigma * x * y)}
#' 
#' \code{GPR <- GPR.polynomial$new(X, y, sigma, p, noise)} with \code{k(x,y) = (x %*% y + sigma)^p}
#'
#' \code{GPR <- GPR.sqrexp$new(X, y, l, noise)} with \code{k(x,y) = exp(-dist(rbind(x, y))^2/(2 * l^2))}
#'
#' \code{GPR <- GPR.gammaexp$new(X, y, gamma, l, noise)} with \code{k(x,y) = exp(-(dist(rbind(x, y)) / l) ^ gamma)}
#'
#' \code{GPR <- GPR.rationalquadratic$new(X, y, alpha, l, noise)} with \code{k(x,y) = (1 + dist(rbind(x, y))^2 / (2 * alpha * l^2))^(-alpha)}
#' 
#' 
#' @importFrom R6 R6Class
#' @name GPR
#' 
#' @examples
#' Hier Beispiele einfügen
#'
#'



GPR <- R6::R6Class("GPR",
                   private = list(
                     .X = NA,
                     .k = NA,
                     .y = NA,
                     .L = NA,
                     .alpha = NA,
                     .noise = NA,
                     .logp = NA
                   ),
                   public = list(
                     initialize = function(X, y, k, noise){
                       stopifnot(is.matrix(X), is.vector(y), is.numeric(y))
                       stopifnot(is.numeric(noise), length(noise) == 1, noise >= 0)
                       stopifnot(is.function(k))
                       private$.X <- X
                       private$.y <- y
                       private$.k <- k
                       private$.noise <- noise
                       n <- ncol(X)
                       K <- outer(1:n, 1:n, function(i,j) k(X[, i, drop = FALSE], X[, j, drop = FALSE]))
                       private$.L <- t(chol(K + noise * diag(n)))
                       private$.alpha <- solve(t(self$L), solve(self$L, y))
                       private$.logp <- -0.5 * self$y %*% self$alpha - sum(log(diag(self$L))) - ncol(self$X) / 2 * log(2 * pi)
                     },
                      predict = function(Xs){
                       
                       #Calculate k* (ks)
                       ks <- sapply(1:ncol(self$X), FUN = function(i) self$k(self$X[, i], Xs))
                       
                       #calculate all other variables directly
                       fs <- ks %*% self$alpha
                       v <- solve(self$L, ks)
                       Vfs <- self$k(Xs, Xs) - v %*% v
                       return(c(fs, Vfs))
                      },
                     posterior_mean = function(testpoints) {
                       outer(1:ncol(testpoints), 1:ncol(self$X), 
                             function(i,j) self$k(testpoints[, i, drop = F], self$X[, j, drop = F])) %*% self$alpha
                     },
                     posterior_covariance = function(testpoints) {
                       len <- ncol(testpoints)
                       cov_test <- outer(1:len, 1:len, 
                              function(i,j) self$k(testpoints[, i, drop = F], testpoints[, j, drop = F]))
                       cov_mixed <- outer(1:len, 1:ncol(self$X), function(i,j) self$k(testpoints[, i, drop = F], self$X[, j, drop = F]))
                       cov_test - cov_mixed %*% sapply(1:len, function(i) solve(t(self$L), solve(self$L, cov_mixed[i, ])))
                     },
                     plot = function(testpoints){
                       dat <- data.frame(x = testpoints, 
                                  y = t(sapply(testpoints, function(x) self$predict(x)[1:2])))
                       ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y.1)) +
                         ggplot2::theme_classic() +
                         ggplot2::scale_y_continuous("output, f(x)") +
                         ggplot2::geom_line() +
                         ggplot2::geom_ribbon(ggplot2::aes(ymin = y.1 - 2*sqrt(pmax(y.2,0)),
                                         ymax = y.1 + 2*sqrt(pmax(y.2,0))), alpha = 0.2) +
                         ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = self$y), 
                                    mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = 4)) +
                         ggplot2::scale_shape_identity()
                     },
                     plot_posterior_draws = function(n, testpoints) {
                       stopifnot(is.matrix(testpoints))
                       len <- ncol(testpoints)
                       mean <- self$posterior_mean(testpoints)
                       covariance <- self$posterior_covariance(testpoints)
                       plot(testpoints, multivariate_normal(len, mean, covariance), type = "l")
                       replicate(n - 1, lines(testpoints, multivariate_normal(len, mean, covariance)))
                     }
                   ),
                   active = list(
                     X = function(value){
                       if(missing(value)){
                         private$.X
                       } else{
                         stop("`$X` is read only", call. = FALSE)
                       }
                     },
                     k = function(value){
                       if(missing(value)){
                         private$.k
                       } else{
                         stop("`$k` is read only", call. = FALSE)
                       }
                     },
                     y = function(value){
                       if(missing(value)){
                         private$.y
                       } else{
                         stop("`$y` is read only", call. = FALSE)
                       }
                     },
                     noise = function(value){
                       if(missing(value)){
                         private$.noise
                       } else{
                         stop("`$noise` is read only", call. = FALSE)
                       }
                     },
                     L = function(value){
                       if(missing(value)){
                         private$.L
                       } else{
                         stop("`$L` is read only", call. = FALSE)
                       }
                     },
                     alpha = function(value){
                       if(missing(value)){
                         private$.alpha
                       } else{
                         stop("`$alpha` is read only", call. = FALSE)
                       }
                     },
                     logp = function(value){
                       if(missing(value)){
                         private$.logp
                       } else{
                         stop("`$logp` is read only", call. = FALSE)
                       }
                     }
                   )
)


GPR.constant <- R6::R6Class("GPR.constant",
                          inherit = GPR,
                          public = list(
                            initialize = function(X, y, c, noise){
                              stopifnot(is.numeric(c), c > 0)
                              k <- function(x, y) constant(x, y, c)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

GPR.linear <- R6::R6Class("GPR.linear", inherit = GPR,
                          public = list(
                            initialize = function(X, y, sigma, noise){
                              stopifnot(length(sigma) == nrow(X))
                              k <- function(x, y) linear(x, y, sigma)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

GPR.polynomial <- R6::R6Class("GPR.polynomial", inherit = GPR,
                              public = list(
                                initialize = function(X, y, sigma, p, noise){
                                  stopifnot(length(sigma) == 1, length(p) == 1)
                                  k <- function(x, y) polynomial(x, y, sigma, p)
                                  super$initialize(X, y, k, noise)
                                }
                              )
)

GPR.sqrexp <-  R6::R6Class("GPR.sqrexp", inherit = GPR,
                           public = list(
                             initialize = function(X, y, l, noise){
                               stopifnot(length(l) == 1)
                               k <- function(x, y) sqrexp(x, y, l)
                               super$initialize(X, y, k, noise)
                             }
                             
                           )
)

GPR.gammaexp <- R6::R6Class("GPR.gammaexp", inherit = GPR,
                          public = list(
                            initialize = function(X, y, gamma, l, noise){
                              stopifnot(length(gamma) == 1, length(l) == 1)
                              k <- function(x, y) gammaexp(x, y, l, gamma)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

GPR.rationalquadratic <- R6::R6Class("GPR.rationalquadratic", inherit = GPR,
                            public = list(
                              initialize = function(X, y, alpha, l, noise){
                                stopifnot(length(alpha) == 1, length(l) == 1)
                                k <- function(x, y) rationalquadratic(x, y, l, alpha)
                                super$initialize(X, y, k, noise)
                              }
                            )
)

# Implementation von Matrixversionen der Kovarianzfunktionen, um Effizienz zu erhöhen
# Bei Eingabe zweier Matrizen gleicher Dimension, soll die Funktion auf die jeweils i-ten Spalten
# für i = 1,...,ncol angewendet werden.
# stopifnots für Dimension einfügen ?!
constant <- function(x, y, c) UseMethod("constant")
constant.matrix <- function(x, y, c) rep(c, ncol(x))
constant.numeric <- function(x, y, c) c

linear <- function(x, y, sigma) UseMethod("linear")
linear.matrix <- function(x, y, sigma) colSums(sigma * x * y)
linear.numeric <- function(x, y, sigma) sum(sigma * x * y)

polynomial <- function(x, y, sigma, p) UseMethod("polynomial")
polynomial.matrix <- function(x, y, sigma, p) (colSums(x * y) + sigma)^p
polynomial.numeric <- function(x, y, sigma, p) (x %*% y + sigma)^p

sqrexp <- function(x, y, l) UseMethod("sqrexp")
sqrexp.matrix <- function(x, y, l) exp(-colSums((x - y)^2)/(2 * l^2))
sqrexp.numeric <- function(x, y, l) exp(-sum((x - y)^2)/(2 * l^2))

gammaexp <- function(x, y, l, gamma) UseMethod("gammaexp")
gammaexp.matrix <- function(x, y, l, gamma) exp(-(sqrt(colSums((x - y)^2))/l)^gamma)
gammaexp.numeric <- function(x, y, l, gamma) exp(-(sqrt(sum((x - y)^2))/l)^gamma)

rationalquadratic <- function(x, y, l, alpha) UseMethod("rationalquadratic")
rationalquadratic.matrix <- function(x, y, l, alpha) (1 + sqrt(colSums((x-y)^2)) / (2 * alpha * l^2))^(-alpha)
rationalquadratic.numeric <- function(x, y, l, alpha) (1 + sqrt(sum((x-y)^2)) / (2 * alpha * l^2))^(-alpha)


cov_dict <- list(
  sqrexp = list(func = function(x, y,l) exp(-sum((x - y)^2)/(2 * l^2)), 
                deriv = function(x, y, l){
                  r <- sqrt(sum((x - y)^2))
                  r^2/l^3*exp(-r^2/(l^2*2))
                }, start = c(1)
          ),
  gammaexp = list(func = function(x, y, gamma, l) exp(-(sqrt(sum((x - y)^2)) / l) ^ gamma), 
                  deriv = function(x, y, gamma, l){
                    r <- sqrt(sum((x - y)^2))
                    c(-exp(-(r/l)^gamma) * (r/l)^gamma * log(r/l), exp(-(r/l)^gamma) * gamma * r^gamma / (l^(gamma + 1)))
                  }, start = c(1, 1)
          ),
  constant = list(func = function(x, y, c) c, deriv = function(x, y, c) 0, start = c(1)
          ),
  linear = list(func = function(x, y, sigma) sum(sigma * x * y), deriv = function(x, y, sigma) sigma, start = c(1)
          ),
  polynomial = list(func = function(x, y, sigma, p) (x %*% y + sigma)^p, 
                  deriv = function(x, y, sigma, p){
                    c(p * (x %*% y + sigma)^(p - 1), (x %*% y + sigma)^p * log((x %*% y + sigma)))
                }, start = c(1, 2)
          )
)

fit <-  function(X, y, noise, cov_names){
  param <- list()
  score <- c()
  for (cov in cov_names){
    usedcov <- cov_dict[[cov]]
    nparam <- length(usedcov$start)
    dens <- function(v){
      n <- ncol(X)
      K <- matrix(0, nrow = n, ncol = n)
      for (i in 1:n) {
        for (j in 1:n) {
          K[i, j] <-  do.call(usedcov$func, as.list(c(X[, i], X[, j], v)))
        }
      }
      
      L <- t(chol(K + noise * diag(n)))
      alpha <- solve(t(L), solve(L, y))
      - 0.5 * y %*% alpha - sum(log(diag(L))) - ncol(X) / 2 * log(2 * pi)
    }
    dens_deriv <- function(v){
      n <- ncol(X)
      K <- matrix(0, nrow = n, ncol = n)
      K_deriv <- array(0, c(n, n, nparam))
      for (i in 1:n) {
        for (j in 1:n) {
          K[i, j] <-  do.call(usedcov$func, as.list(c(X[, i], X[, j], v)))
          K_deriv[i, j,] <- do.call(usedcov$deriv, as.list(c(X[, i], X[, j], v)))
        }
      }
      
      K_inv <- solve(K)
      alpha <- K_inv %*% y
      vapply(1:nparam, function(i) 0.5 * sum(diag(alpha %*% t(alpha) - K_inv) %*% K_deriv[,,i]), numeric(1))
    }
    #Switch optim method if parameter is one dimensional
    if (nparam == 1) p <- optim(usedcov$start, dens, gr = dens_deriv, method = "Brent", 
                                upper = 10, lower = -10, control = list(fnscale = -1))
    else p <- optim(usedcov$start, dens, gr = dens_deriv, method = "BFGS", 
                   control = list(fnscale = -1))
    
    param <- append(param, list(p$par))
    score <- c(score, p$value)
  }
  return(list(par = param[[which.max(score)]], cov = cov_names[[which.max(score)]], score = score))
}


X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
noise <- 0.5
y <- c(0.1*X^3 + rnorm(length(X),0, 1))

Gaussian <- GPR.constant$new(X, y, 1, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian <- GPR.linear$new(X, y, 1, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian <- GPR.polynomial$new(X, y, 1, 3, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian <- GPR.rationalquadratic$new(X, y, 1, 1.5, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian <- GPR.sqrexp$new(X, y, 1, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian <- GPR.gammaexp$new(X, y, 1, 1.5, noise)
Gaussian$plot(seq(-5,5, by = 0.1))


z <- fit(X,y,noise,list("sqrexp", "gammaexp"))

print(z)
Gaussian <- GPR$new(X, y, function(x,y) do.call(cov_dict[[z$cov]]$func, append(list(x,y),z$par)), noise)
#Gaussian <- GPR$new(X, y, function(a,b) cov_dict[[z$cov]]$func(a,b,z$par), noise)
#Gaussian$plot(seq(-5,5, by = 0.1))


X <- matrix(seq(-5,5,by = 0.5), nrow = 1)
noise <- 0.5
y <- c(0.1*X^3 + rnorm(length(X), 0, 1))
Gaussian <- GPR.gammaexp$new(X, y, 1, 1.5, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian$plot_posterior_draws(10, matrix(seq(-5,5, by = 0.1), nrow = 1))
