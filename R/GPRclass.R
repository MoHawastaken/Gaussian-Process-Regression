#'  Predictions and Plots for Gaussian Process Regression
#'
#'  Implements gaussian processes and gives tools for gaussian process regression and classification problems for given testpoints including clear plots of the results.
#' 
#' @usage 
#' \preformatted{GPR <- GPR$new(X, y, cov_fun, noise)
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
#'   \code{cov_fun} the chosen covariance function of the gaussian process
#' 
#'   \code{noise} the inflicted noise of the observations
#' 
#'   \code{X*} a numeric vector as the test input
#' 
#'   \code{testpoints} a matrix of testpoints
#'   
#' @section Methods:
#' \code{$predict()} returns a numeric vector of the expected value of the underlying function f and its variance for the test input
#' 
#' \code{$plot()} displays the results of the \code{predict} function for all testpoints and confidence regions of two standard deviations in a transparent plot
#' 
#'
#' @section Subclasses:
#'
#' GPR has several subclasses where a covariance function k(x,y) is given. The following subclasses are implemented:
#' 
#' \code{GPR <- GPR.constant$new(X, y, c, noise)} with \code{k(x,y) = c}
#' 
#' \code{GPR <- GPR.linear$new(X, y, cov_Fun, noise)} with \code{k(x,y) = sum(sigma * x * y)}
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
#' @importFrom R6 R6Class
#' @name GPR
#' @section Details:
#' If own covariance functions are used with GPR, they need to be vectorized.
#' 
#' @examples
#' Hier Beispiele einfuegen
#'
#' @references Rasmussen, Carl E.; Williams, Christopher K. I. (2006).	Gaussian processes for machine learning
NULL

#' @export
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
                 #K <- outer(1:n, 1:n, function(i,j) k(X[, i, drop = FALSE], X[, j, drop = FALSE]))
                 K <- covariance_matrix(X, X, k)
                 #Pruefe, ob alle Hauptminoren positiv sind.
                 if (min(sapply(1:n, function(i) det((K + noise * diag(n))[1:i, 1:i, drop = F]))) <= 0) {
                   stop("Inputs lead to non positive definite covariance matrix. Try using a larger noise or a smaller lengthscale.")
                 }
                 private$.L <- t(chol(K + noise * diag(n)))
                 private$.alpha <- solve(t(self$L), solve(self$L, y))
                 private$.logp <- -0.5 * self$y %*% self$alpha - sum(log(diag(self$L))) - ncol(self$X) / 2 * log(2 * pi)
               },
               predict = function(X_star, pointwise_var = FALSE){
                 stopifnot(is.numeric(X_star), length(X_star) %% nrow(self$X) == 0)
                 if (is.null(dim(X_star))) {
                   dim(X_star) <- c(nrow(self$X), length(X_star)/nrow(self$X))
                 }
                 K_star <- covariance_matrix(self$X, X_star, self$k)
                 posterior_mean <- t(K_star) %*% self$alpha
                 v <- solve(self$L, K_star)
                 if (pointwise_var) {
                   posterior_variance <- self$k(X_star, X_star) - colSums(v * v)
                   return(cbind(posterior_mean, posterior_variance))
                 } else {
                   posterior_variance <- covariance_matrix(X_star, X_star, self$k) - t(v) %*% v
                   return(list(posterior_mean, posterior_variance))
                 }
               },
               plot = function(testpoints){
                 y <- unname(self$predict(testpoints, pointwise_var = TRUE))
                 dat <- data.frame(x = testpoints, y = y)
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
                 predictions <- self$predict(testpoints)
                 y <- cbind(predictions[[1]], diag(predictions[[2]]))
                 z <- multivariate_normal(n, predictions[[1]], predictions[[2]])
                 dat <- data.frame(x = testpoints, y = y, z = z)
                 dat1 <- tidyr::gather(dat, -(x:y.2), key = "variable", value = "value")
                 ggplot2::ggplot(dat1, ggplot2::aes(x = x, y = value, colour = variable)) +
                   ggplot2::theme_classic() +
                   ggplot2::scale_y_continuous("Random functions drawn from posterior") +
                   ggplot2::geom_line() +
                   ggplot2::geom_ribbon(inherit.aes = F, mapping = ggplot2::aes(x = testpoints, ymin = y.1 - 2*sqrt(pmax(y.2,0)),
                                                     ymax = y.1 + 2*sqrt(pmax(y.2,0))), data = dat, alpha = 0.3) 
                   #ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = self$y), 
                    #              mapping = ggplot2::aes(x = xpoints, y = ypoints))
              },
              plot_posterior_variance = function(where, limits = c(min(self$X), max(self$X)), subdivisions = 100L) {
                x <- seq(limits[1], limits[2], length.out = subdivisions)
                len <- length(where)
                y <- self$predict(c(where, x))[[2]][(len + 1):(len + subdivisions), 1:len]
                dat <- data.frame(x = x, y = y)
                names(dat) <- c("x", as.character(where))
                dat <- tidyr::gather(dat, -x, key = "z", value = "value")
                ggplot2::ggplot(dat, ggplot2::aes(x = x, y = value, colour = z)) +
                  ggplot2::theme_classic() +
                  ggplot2::scale_y_continuous("Posterior covariance of f(x) with f(z)") +
                  ggplot2::geom_line()
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

#' @export
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

#' @export
GPR.linear <- R6::R6Class("GPR.linear", inherit = GPR,
                          public = list(
                            initialize = function(X, y, sigma, noise){
                              stopifnot(length(sigma) == nrow(X))
                              k <- function(x, y) linear(x, y, sigma)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

#' @export
GPR.polynomial <- R6::R6Class("GPR.polynomial", inherit = GPR,
                              public = list(
                                initialize = function(X, y, sigma, p, noise){
                                  stopifnot(length(sigma) == 1, length(p) == 1)
                                  k <- function(x, y) polynomial(x, y, sigma, p)
                                  super$initialize(X, y, k, noise)
                                }
                              )
)

#' @export
GPR.sqrexp <-  R6::R6Class("GPR.sqrexp", inherit = GPR,
                           public = list(
                             initialize = function(X, y, l, noise){
                               stopifnot(length(l) == 1)
                               k <- function(x, y) sqrexp(x, y, l)
                               super$initialize(X, y, k, noise)
                             }
                             
                           )
)

#' @export
GPR.gammaexp <- R6::R6Class("GPR.gammaexp", inherit = GPR,
                          public = list(
                            initialize = function(X, y, gamma, l, noise){
                              stopifnot(length(gamma) == 1, length(l) == 1)
                              k <- function(x, y) gammaexp(x, y, l, gamma)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

#' @export
GPR.rationalquadratic <- R6::R6Class("GPR.rationalquadratic", inherit = GPR,
                            public = list(
                              initialize = function(X, y, alpha, l, noise){
                                stopifnot(length(alpha) == 1, length(l) == 1)
                                k <- function(x, y) rationalquadratic(x, y, l, alpha)
                                super$initialize(X, y, k, noise)
                              }
                            )
)

# Funktion, um Kovarianzmatrix zu berechnen. Die verwendete Kovarianzfunktion muss bei Eingabe zweier Matrizen 
# gleicher Dimension die Werte bei Anwendung auf die jeweils i-ten Spalten für i = 1,...,ncol zurueckgeben.
covariance_matrix <- function(A, B, covariance_function) {
  outer(1:ncol(A), 1:ncol(B), function(i, j) covariance_function(A[, i, drop = F], B[, j, drop = F]))
}
multivariate_normal <- function(n, mean, covariance) {
  stopifnot(is.numeric(mean), is.numeric(covariance), length(mean) == nrow(covariance))
  degenerate <- diag(covariance) < 0.05
  L <- t(chol(covariance[!degenerate, !degenerate]))
  out <- matrix(0, nrow = nrow(covariance), ncol = n)
  out[!degenerate] <- L%*%matrix(rnorm(n*nrow(L), 0, 1), nrow = nrow(L))
  out + c(mean) 
}

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
rationalquadratic.matrix <- function(x, y, l, alpha) (1 + sqrt(colSums((x - y)^2)) / (2 * alpha * l^2))^(-alpha)
rationalquadratic.numeric <- function(x, y, l, alpha) (1 + sqrt(sum((x - y)^2)) / (2 * alpha * l^2))^(-alpha)

X <- matrix(seq(-5,5,by = 0.5), nrow = 1)
noise <- 0.5
y <- c(0.1*X^3 + rnorm(length(X), 0, 1))
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

Gaussian <- GPR.gammaexp$new(X, y, 1, 1.5, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
Gaussian$plot_posterior_draws(10, seq(-5,5, by = 0.1))

X <- matrix(seq(-5,5, by = 3), nrow = 1)
noise <- 0
y <- c(0.1*X^3 + rnorm(length(X), 0, 1))
Gaussian <- GPR.sqrexp$new(X, y, 0.5, 0)
Gaussian$plot(seq(-6,6, by = 0.1))
Gaussian$plot_posterior_draws(3, seq(-6,6, by = 0.2))