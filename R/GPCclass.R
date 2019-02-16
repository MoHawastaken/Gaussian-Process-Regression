#'  Predictions and Plots for Gaussian process classification
#'
#'  Implements a gaussian process and gives tools to predict and plot its values for given testpoints
#' 
#'
#' @usage \preformatted{GPC <- GPC$new(X, y, cov_Fun, noise)
#'
#'
#' GPC$predict(X*)
#' GPC$plot(testpoints)
#'}
#' @section Arguments:
#' 
#'   \code{X} matrix of inputs
#'
#'   \code{y} numeric vector of targets
#' 
#'   \code{cov_fun} the chosen covariance function of the gaussian process
#' 
#'   \code{noise} the inflicted noise of the observations
#' 
#'   \code{X*} a numeric vector as the test input
#' 
#'   \code{testpoints} a matrix of testpoints
#'   
#' @section Methods:
#' 
#' \code{$predict()} returns a numeric vector of the expected value of the underlying 
#' function f and their variance for the test input
#' 
#' \code{$plot()} displays the results of the predict function for all testpoints in a nice plot
#'
#' 
#' 
#' @examples
#' Hier Beispiele einfügen
#'
#'
#' @importFrom R6 R6Class
#' @name GPC
NULL

#' @export
GPC <- R6::R6Class("GPC",
                   private = list(
                     .X = NA,
                     .k = NA,
                     .y = NA,
                     .f_hat = NA,
                     .L = NA,
                     .logq = NA,
                     .sigmoid = function(x) 1/(1 + exp(-x)) #used sigmoid function
                   ),
                   public = list(
                     initialize = function(X, y, k, epsilon){
                       stopifnot(is.matrix(X), is.vector(y), is.numeric(y), length(y) == ncol(X))
                       stopifnot(is.numeric(epsilon), epsilon > 0, is.function(k))
                       n <- length(y)
                       K <- covariance_matrix(X,X,k)
                       f <- rep(0, n)
                       it <- 0
                       while (TRUE) {
                         it <- it + 1
                         P <- private$.sigmoid(f)
                         W <- (1 - P) * P
                         L <- t(chol(diag(n) + (sqrt(W) %o% sqrt(W)) * K))
                         b <- W * f + (y + 1)/2 - P
                         intermediate <- solve(L, sqrt(W) * (K %*% b))
                         intermediate <- solve(t(L), intermediate)
                         a <- b - sqrt(W) * intermediate
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
                       P <- private$.sigmoid(f)
                       W <- (1 - P) * P
                       private$.f_hat <- f
                       private$.L <- t(chol(diag(n) + (sqrt(W) %o% sqrt(W)) * K))
                       private$.logq <- objective - sum(diag(self$L))
                       private$.X <- X
                       private$.y <- y
                       private$.k <- k
                     },
                     predict_class = function(Xs){
                       P <- private$.sigmoid(self$f_hat)
                       W <- P * (1 - P)
                       ks <- sapply(1:ncol(self$X), FUN = function(i) self$k(self$X[, i], Xs))
                       fs_bar <- sum(ks * ((self$y + 1)/2 - P))
                       v <- solve(self$L, (sqrt(W) * ks))
                       Vfs <- self$k(Xs, Xs) - sum(v * v)
                       hilfs_func <- function(z) private$.sigmoid(z) * dnorm(z, mean = fs_bar, sd = Vfs)
                       P_hat <- integrate(hilfs_func, -Inf, Inf)[1]
                       return(P_hat$value)
                     },
                     plot = function(testpoints){
                       dat <- data.frame(x = testpoints, 
                                         y = sapply(testpoints, function(x) self$predict_class(x)))
                       ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
                         ggplot2::theme_classic() +
                         ggplot2::scale_y_continuous("output, p(y = 1| x)") +
                         ggplot2::geom_line() +
                         #ggplot2::geom_ribbon(ggplot2::aes(ymin = y.1 - 2*sqrt(pmax(y.2,0)),
                         #                                  ymax = y.1 + 2*sqrt(pmax(y.2,0))), alpha = 0.2) +
                         ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = pmax(0,self$y)), 
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

#section for testing:

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
