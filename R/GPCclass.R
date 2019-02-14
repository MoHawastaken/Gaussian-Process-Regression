#'  Predictions and Plots for Gaussian process regression
#'
#'  Implements a gaussian process and gives tools to predict and plot its values for given testpoints
#' 
#'
#' @section Usage: 
#' 
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
#' @section Methods:
#' 
#' \code{$predict()} returns a numeric vector of the expected value of the underlying function f and their variance for the test input
#' 
#' \code{$plot()} displays the results of the predict function for all testpoints in a nice plot
#' 
#'
#' @section Subclasses:
#' 
#' GPR has several subclasses where a covarianz function k(x,y) is given. The following subclasses are implemented
#' 
#' \code{GPR <- GPR.constant$new(X, y, c, noise)} with \code{k(x,y) = c}
#' 
#' \code{GPR <- GPR.linear$new(X, y, sigma, noise)} with \code{k(x,y) = sum(sigma * x * y)}
#' 
#' \code{GPR <- GPR.polynomial$new(X, y, sigma, p, noise)} with \code{k(x,y) = (x \%*\% y + sigma)^p}
#'
#' \code{GPR <- GPR.sqrexp$new(X, y, l, noise)} with \code{k(x,y) = exp(-dist(rbind(x, y))^2/(2 * l^2))}
#'
#' \code{GPR <- GPR.gammaexp$new(X, y, gamma, l, noise)} with \code{k(x,y) = exp(-(dist(rbind(x, y)) / l) ^ gamma)}
#'
#' \code{GPR <- GPR.rationalquadratic$new(X, y, alpha, l, noise)} with \code{k(x,y) = (1 + dist(rbind(x, y))^2 / (2 * alpha * l^2))^(-alpha)}
#' 
#'
#' 
#' 
#' @examples
#' Hier Beispiele einfügen
#'
#'
#' @importFrom R6 R6Class
#' @name GPR


GPC <- R6::R6Class("GPR",
                   private = list(
                     .X = NA,
                     .k = NA,
                     .y = NA,
                     .f_hat = NA,
                     .L = NA,
                     .logq = NA
                   ),
                   public = list(
                     initialize = function(X, y, k, epsilon){
                       stopifnot(is.matrix(X), is.vector(y), is.numeric(y), length(y) == ncol(X))
                       stopifnot(is.numeric(epsilon), epsilon > 0, is.function(k))
                       n <- length(y)
                       K <- matrix(0, nrow = n, ncol = n)
                       for (i in 1:n) {
                         for (j in 1:n) {
                           K[i, j] <-  k(X[, i], X[, j])
                         }
                       }
                       f <- rep(0, n)
                       it <- 0
                       while (TRUE) {
                         it <- it + 1
                         Pi <- 1/(1 + exp(-f))
                         W <- (1 - Pi)*Pi
                         L <- t(chol(diag(n) + (sqrt(W) %o% sqrt(W)) * K))
                         b <- W*f + (y + 1)/2 - Pi
                         intermediate <- solve(L, sqrt(W)*(K %*% b))
                         intermediate <- solve(t(L), intermediate)
                         a <- b - sqrt(W)*intermediate
                         f <- c(K %*% a)
                         objective <- -sum(a * f)/2 - sum(log(1 + exp(-y * f)))
                         if (it > 1) {
                           if (abs(objective - last_objective) < epsilon) {
                             break
                           } else if (least_objective + 10 < objective) {
                             stop("Apparently does not converge.")
                           }
                         } else {
                           least_objective <- objective
                         }
                         last_objective <- objective
                       }
                       print(sprintf("Convergence after %s iterations", it))
                       Pi <- 1/(1 + exp(-f))
                       W <- (1 - Pi)*Pi
                       private$.f_hat <- f
                       private$.L <- t(chol(diag(n) + (sqrt(W) %o% sqrt(W)) * K))
                       private$.logq <- objective - sum(diag(self$L))
                       private$.X <- X
                       private$.y <- y
                       private$.k <- k
                     },
                     predict_class = function(Xs){
                       Pi <- 1/(1 + exp(-self$f_hat))
                       W <- Pi*(1 - Pi)
                       ks <- sapply(1:ncol(self$X), FUN = function(i) self$k(self$X[, i], Xs))
                       fs_bar <- sum(ks * ((self$y + 1)/2 - Pi))
                       v <- solve(self$L, (sqrt(W) * ks))
                       Vfs <- self$k(Xs, Xs) - sum(v * v)
                       hilfs_func <- function(z) 1/(1 + exp(-z))  * (1/sqrt(2*pi*Vfs)) * exp(-(z-fs_bar)^2/(2*Vfs))
                       PIs_hat <- integrate(hilfs_func, -Inf, Inf)[1]
                       return(PIs_hat$value)
                     },
                     plot = function(testpoints){
                       dat <- data.frame(x = testpoints, 
                                         y = sapply(testpoints, function(x) self$predict_class(x)))
                       ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
                         ggplot2::theme_classic() +
                         ggplot2::scale_y_continuous("output, p(y = 1| x)") +
                         ggplot2::geom_line() +
                         #ggplot2::geom_ribbon(ggplot2::aes(ymin = y.1 - 2*sqrt(max(y.2,0)),
                         #                                  ymax = y.1 + 2*sqrt(max(y.2,0))), alpha = 0.2) +
                         ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = self$y), 
                                             mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = 4)) +
                         ggplot2::scale_shape_identity()
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
                     f_hat = function(value){
                       if(missing(value)){
                         private$.f_hat
                       } else{
                         stop("`$f_hat` is read only", call. = FALSE)
                       }
                     },
                     L = function(value){
                       if(missing(value)){
                         private$.L
                       } else{
                         stop("`$L` is read only", call. = FALSE)
                       }
                     },
                     logq = function(value){
                       if(missing(value)){
                         private$.logq
                       } else{
                         stop("`$logq` is read only", call. = FALSE)
                       }
                     }
                   )
)

X <- matrix(seq(-1,1,by = 0.1), nrow = 1)
y <- 2*as.integer(X > 0) - 1
kappa <- function(x,y) exp(-3*(x - y)^2)
gaussian_classifier <- GPC$new(X, y, kappa, 1e-5)
gaussian_classifier$plot(seq(-2,2, by = 0.1))

X <- matrix(c(seq(-1, -0.1, by = 0.1), seq(0, 1, by = 0.2)), nrow = 1)
y <- 2*as.integer(X > 0) - 1
kappa <- function(x,y) exp(-3*(x - y)^2)
gaussian_classifier <- GPC$new(X, y, kappa, 1e-5)
gaussian_classifier$plot(seq(-2,2, by = 0.1))
