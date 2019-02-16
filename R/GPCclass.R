#'  Predictions and Plots for Gaussian process classification
#'
#'  Implements a gaussian process for classification and gives tools 
#'  to predict and plot its values for given testpoints
#' 
#'
#' @usage \preformatted{GPC <- GPC$new(X, y, cov_fun, noise)
#'
#' GPC$predict(X*)
#' GPC$plot(testpoints)
#'}
#' @section Arguments:
#' 
#'   \code{X} matrix of inputs
#'
#'   \code{y} numeric vector of targets with values in \{-1,1\}
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
#' function f
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
#' @references Rasmussen, Carl E.; Williams, Christopher K. I. (2006).	Gaussian processes for machine learning
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
                       f <- rep(0, n) #starting value
                       it <- 0 #number of iterations
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
                       message(sprintf("Convergence after %s iterations", it))
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
                     predict_class2 = function(X_star) {
                       P <- private$.sigmoid(self$f_hat)
                       W <- P * (1 - P)
                       K_star <- covariance_matrix(self$X, X_star, self$k)
                       fs_bar <- t(K_star) %*% ((self$y + 1)/2 - P)
                       v <- solve(self$L, (sqrt(W) * K_star))
                       Vfs <- self$k(X_star, X_star) - colSums(v * v)
                       sapply(1:length(Vfs), function(i) integrate(function(z) 
                                private$.sigmoid(z) * dnorm(z, mean = fs_bar[i, 1], sd = Vfs[i]), -Inf, Inf)$value)
                     },
                     plot = function(testpoints){
                       if(is.vector(testpoints)){
                         dat <- data.frame(x = testpoints, 
                                           y = sapply(testpoints, function(x) self$predict_class(x)))
                         ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
                           ggplot2::theme_classic() +
                           ggplot2::scale_y_continuous("output, p(y = 1| x)") +
                           ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = factor(as.integer(y > 0.5))), 
                                              height = Inf, alpha = 0.5) +
                           ggplot2::guides(fill = ggplot2::guide_legend(title = "Labels")) +
                           ggplot2::scale_fill_manual(values = c("red", "blue")) +
                           ggplot2::geom_line() +
                           #ggplot2::geom_ribbon(ggplot2::aes(ymin = y.1 - 2*sqrt(pmax(y.2,0)),
                           #                                  ymax = y.1 + 2*sqrt(pmax(y.2,0))), alpha = 0.2) +
                           ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = pmax(0,self$y)), 
                                               mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = 4)) +
                           ggplot2::scale_shape_identity() 
                       }
                       else if(nrow(testpoints) == 2){
                         #dat <- data.frame(x.1 = testpoints[1,], x.2 = testpoints[2,],
                            #y = apply(testpoints, 2, function(x) {
                             # 2*as.integer(self$predict_class(x) >= 0.5) - 1}))
                        dat <- data.frame(x.1 = testpoints[1,], x.2 = testpoints[2,],
                                  y = 2*as.integer(self$predict_class2(testpoints) >= 0.5) - 1)
                         ggplot2::ggplot(dat, inherit.aes = F, ggplot2::aes(x = x.1, y = x.2, fill = factor(y))) +
                           ggplot2::theme_classic() +
                           ggplot2::scale_y_continuous(expression("x_2")) +
                           ggplot2::scale_x_continuous(expression("x_1")) +
                           ggplot2::geom_tile() + 
                           ggplot2::scale_fill_manual(values = c("red", "blue")) +
                           ggplot2::guides(fill = ggplot2::guide_legend(title = "Labels")) +
                           ggplot2::geom_point(inherit.aes = F, data = data.frame(xpoints = c(self$X[1,]), ypoints = c(self$X[2,])), 
                                               mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = factor(self$y))) +
                           ggplot2::scale_shape_manual(values = c(4, 2)) +
                           ggplot2::guides(shape = ggplot2::guide_legend(title = "Testpoints")) +
                           ggplot2::scale_color_manual(values = c("red", "blue"))
                         }
                       else warning("Plot function not available for this dimension")
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

s <- seq(-1, 1, by = 0.5)
X <- matrix(c(rep(s, each = length(s)), rep(s, times = length(s))), nrow = 2, byrow = T)
y <- 2*as.integer(X[1, ] > X[2, ]) - 1
kappa <- function(x,y) sqrexp(x,y,l=1)
gaussian_classifier <- GPC$new(X, y, kappa, 1e-5)
s <- seq(-1, 1, by = 0.1)
testpoints <- matrix(c(rep(s, each = length(s)), rep(s, times = length(s))), nrow = 2, byrow = T)
gaussian_classifier$plot(testpoints)

n <- 10
X <- cbind(multivariate_normal(n, c(0.5,0.5), diag(c(0.1,0.1))), multivariate_normal(n, c(-0.5,-0.5), diag(c(0.1,0.1))))
plot(X[1,], X[2, ])
y <- rep(c(1,-1), each = n)
kappa <- function(x,y) sqrexp(x,y,l=1)
gaussian_classifier <- GPC$new(X, y, kappa, 1e-5)
s <- seq(-1, 1, by = 0.1)
testpoints <- matrix(c(rep(s, each = length(s)), rep(s, times = length(s))), nrow = 2, byrow = T)
gaussian_classifier$plot(testpoints)


