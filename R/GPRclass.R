#'  Predictions and Plots for Gaussian Process Regression
#'
#'  Implements gaussian processes and gives tools for gaussian process
#'  regression problems for given testpoints including clear plots of the
#'  results and optimization of hyperparameters.
#' 
#' @usage 
#' \preformatted{gpr_object <- GPR$new(X, y, noise = 0, k = fit(X, y, noise, cov_names)$func, cov_names = names(cov_dict))
#'
#'
#' gpr_object$predict(X_star, pointwise_var = TRUE)
#' gpr_object$plot(limits, length.out = 200L)
#' gpr_object$plot_posterior_draws(n = 5, limits, length.out = 200L)
#' gpr_object$plot_posterior_variance(where, limits, length.out = 200L)
#'}

#' @section Arguments:
#' 
#'   \code{X}          matrix of inputs
#'
#'   \code{y}          numeric vector of targets
#'   
#'   \code{noise}      the inflicted noise of the observations
#' 
#'   \code{cov_fun}    the chosen covariance function of the gaussian process (optional)
#'   
#'   \code{cov_names}  a list of given covariance functions to optimize over; relevant if \code{cov_fun} is not given (optional)
#' 
#'   \code{X_star}     a matrix of test inputs, where each column is interpreted as one observation
#' 
#'   \code{limits}     a vector of length 2, giving the lower and upper bound of the plot; the default is the extended range of X
#'   
#'   \code{length.out} the number of subdivisions for the plots
#'   
#' @section Methods:
#' \code{$predict(X_star, pointwise_var = TRUE)}  
#' returns a matrix of the expected value of the underlying function f and its variance for the
#' testpoints. If the input is a vector of length n, predict will interpret it
#' as n testpoints. Otherwise each column of the input matrix is interpreted as
#' a testpoint. For \code{pointwise_var = FALSE} , predict will return the
#' predicted covariance matrix cov(X_star, X_star) instead of only its diagonal.
#' 
#' \code{$plot(limits, length.out = 200L)}
#' displays the results of the \code{predict} function on the interval given by \code{limits} and confidence
#' regions of two standard deviations in a transparent plot
#' 
#' \code{$plot_posterior_draws(n = 5, limits, length.out = 200L)}
#' plots n random functions drawn from the posterior distribution of the underlying Gaussian
#' process on the interval given by \code{limits}
#' 
#' \code{$plot_posterior_variance(where, limits, length.out = 200L)}
#' visualizes the posterior covariance of the Gaussian process with specified points
#' \code{where}
#'
#' @section Subclasses:
#'
#' GPR has several subclasses where a covariance function k(x,y) is given. The following subclasses are implemented:
#' 
#' \code{GPR <- GPR.constant$new(X, y, noise, c)}                 with \code{k(x,y) = c}
#' 
#' \code{GPR <- GPR.linear$new(X, y, noise, sigma)}               with \code{k(x,y) = sum(sigma * x * y)}
#' 
#' \code{GPR <- GPR.polynomial$new(X, y, noise, sigma, p)}        with \code{k(x,y) = (x \%*\% y + sigma)^p}
#'
#' \code{GPR <- GPR.sqrexp$new(X, y, noise, l)}                   with \code{k(x,y) = exp(-dist(rbind(x, y))^2/(2 * l^2))}
#'
#' \code{GPR <- GPR.gammaexp$new(X, y, noise, gamma, l)}          with \code{k(x,y) = exp(-(dist(rbind(x, y)) / l) ^ gamma)}
#'
#' \code{GPR <- GPR.rationalquadratic$new(X, y, noise, alpha, l)} with \code{k(x,y) = (1 + dist(rbind(x, y))^2 / (2 * alpha * l^2))^(-alpha)}
#' 
#'
#' 
#' @importFrom R6 R6Class
#' @name GPR
#' @section Details:
#' If own covariance functions are used with GPR, they need to be vectorized.
#' The plot functions are only implemented for one dimensional data. Although
#' this class allows customary use with its amount of optional parameters, since
#' the default values are set by the optimal choice of hyperparameters and
#' preimplemented covariance functions, the class is designed for well fitted
#' and easy use by only giving the necessary inputs \code{X, y, noise}.
#' 
#' @examples
#' X <- matrix(seq(-5,5,by = 1), nrow = 1)
#' y <- c(0.1*X^3 + rnorm(length(X), 0, 1))
#' 
#' #using optimal parameters for gammaexp
#' Gaussian <- GPR.gammaexp$new(X, y, noise = 0.1)
#' # == GPR$new(X, y, noise = 0.1, cov_names = "gammaexp")
#' Gaussian$plot()
#' Gaussian$plot_posterior_draws()
#' Gaussian$plot_posterior_variance(seq(-5, 5, by = 3))
#' 
#' #optimizing hyperparameters over all covariance functions
#' Gaussian <- GPR$new(X, y, noise = 0.1)
#' Gaussian$plot()
#' Gaussian$plot_posterior_draws()
#' Gaussian$plot_posterior_variance(seq(-5, 5, by = 3))
#' 
#' X <- matrix(seq(-5,5, by = 1), nrow = 1)
#' noise <- 0
#' y <- c(0.1*X^3 + rnorm(length(X), 0, 1))
#' # explicit use of squared exponential with l = 0.5
#' Gaussian <- GPR.sqrexp$new(X, y, l = 0.5, noise = 0.5)
#' Gaussian$plot()
#' Gaussian$plot_posterior_draws(3)
#' Gaussian$plot_posterior_variance(seq(-5, 5, by = 3))
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
               initialize = function(X, y, noise = 0, k = fit(X, y, noise, cov_names)$func,
                                     cov_names = names(cov_dict)){
                 stopifnot(is.numeric(X), is.vector(y), is.numeric(y))
                 stopifnot(is.numeric(noise), length(noise) == 1, noise >= 0)
                 # Ist Input X ein Vektor, wird dieser als einzeilige Matrix behandelt.
                 if (!is.matrix(X)) dim(X) <- c(1, length(X))
                 stopifnot(length(y) == ncol(X), is.function(k))
                 private$.X <- X
                 private$.y <- y
                 private$.k <- k
                 private$.noise <- noise
                 n <- ncol(X)
                 K <- covariance_matrix(X, X, k)
                 #Pruefe, ob alle Hauptminoren positiv sind.
                 if (min(sapply(1:n, function(i) det((K + noise * diag(n))[1:i, 1:i, drop = F]))) <= 0) {
                   stop("Inputs lead to non positive definite covariance matrix. Try using a larger noise or a smaller lengthscale.")
                 }
                 if(class(try(solve(K + noise * diag(n)),silent=T)) != "matrix"){
                   stop("K(X,X) + noise * I is not invertible, the algorithm is not defined for this case.")
                 }
                 private$.L <- t(chol(K + noise * diag(n)))
                 private$.alpha <- solve(t(self$L), solve(self$L, y))
                 private$.logp <- -0.5 * self$y %*% self$alpha - sum(log(diag(self$L))) - ncol(self$X) / 2 * log(2 * pi)
               },
               predict = function(X_star, pointwise_var = TRUE){
                 stopifnot(is.numeric(X_star), length(X_star) %% nrow(self$X) == 0)
                 if (is.null(dim(X_star))) {
                   dim(X_star) <- c(nrow(self$X), length(X_star)/nrow(self$X))
                 }
                 K_star <- covariance_matrix(self$X, X_star, self$k)
                 posterior_mean <- t(K_star) %*% self$alpha
                 v <- solve(self$L, K_star)
                 if (pointwise_var) {
                   posterior_variance <- self$k(X_star, X_star) - colSums(v * v)
                   return(cbind(posterior_mean, posterior_variance, deparse.level = 0))
                 } else {
                   posterior_variance <- covariance_matrix(X_star, X_star, self$k) - t(v) %*% v
                   return(list(posterior_mean, posterior_variance))
                 }
               },
               plot = function(limits = c(expand_range(self$X)[1], expand_range(self$X)[2]), length.out = 200L){
                 if (nrow(self$X) > 1) {
                   message("No plot method for multidimensional data.")
                   return
                 }
                 testpoints <- seq(limits[[1]], limits[[2]], length.out = length.out)
                 y <- self$predict(testpoints, pointwise_var = TRUE)
                 dat <- data.frame(x = testpoints, y = y)
                 g <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y.1)) +
                   ggplot2::theme_classic() +
                   ggplot2::scale_y_continuous("output, f(x)") +
                   ggplot2::geom_line() +
                   ggplot2::geom_ribbon(ggplot2::aes(ymin = y.1 - 2*sqrt(pmax(y.2,0)),
                                   ymax = y.1 + 2 * sqrt(pmax(y.2,0))), alpha = 0.2) +
                   ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = self$y), 
                              mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = 4)) +
                   ggplot2::scale_shape_identity()
                  list(plot = g, pred = y)
               },
               plot_posterior_draws = function(n = 5, limits = c(expand_range(self$X)[1], expand_range(self$X)[2]),
                                                                       length.out = 200L) {
                 if (nrow(self$X) > 1) {
                   message("No plot method for multidimensional data.")
                   return
                 }
                 testpoints <- seq(limits[[1]], limits[[2]], length.out = length.out)
                 predictions <- self$predict(testpoints, pointwise_var = FALSE)
                 y <- cbind(predictions[[1]], diag(predictions[[2]]))
                 z <- multivariate_normal(n, predictions[[1]], predictions[[2]])
                 dat <- data.frame(x = testpoints, y = y, z = z)
                 dat1 <- tidyr::gather(dat, -(x:y.2), key = "variable", value = "value")
                 ggplot2::ggplot(dat1, ggplot2::aes(x = x, y = value, colour = variable)) +
                   ggplot2::theme_classic() +
                   ggplot2::scale_y_continuous("Random functions drawn from posterior") +
                   ggplot2::geom_line() +
                   ggplot2::geom_ribbon(inherit.aes = F, data = dat, mapping = ggplot2::aes(x = testpoints, 
                                    ymin = y.1 - 2*sqrt(pmax(y.2,0)), ymax = y.1 + 2*sqrt(pmax(y.2,0))), alpha = 0.3) +
                   ggplot2::guides(colour = FALSE)
                   #ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = self$y), 
                  #               mapping = ggplot2::aes(x = xpoints, y = ypoints))
              },
              plot_posterior_variance = function(where, limits = c(expand_range(self$X)[1], expand_range(self$X)[2]),
                                                                         length.out = 200L) {
                if (nrow(self$X) > 1) {
                  message("No plot method for multidimensional data.")
                  return
                }
                testpoints <- seq(limits[[1]], limits[[2]], length.out = length.out)
                len <- length(where)
                y <- self$predict(c(where, testpoints), pointwise_var = FALSE)[[2]][(len + 1):(len + length(testpoints)), 1:len]
                dat <- data.frame(x = testpoints, y = y)
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
                            initialize = function(X, y, noise, c = fit(X,y,noise,"constant")$par){
                              stopifnot(is.numeric(c), c > 0)
                              k <- function(x, y) constant(x, y, c)
                              super$initialize(X, y, noise, k)
                            }
                          )
)

#' @export
GPR.linear <- R6::R6Class("GPR.linear", inherit = GPR,
                          public = list(
                            initialize = function(X, y, noise, sigma = fit(X,y,noise,"linear")$par){
                              stopifnot(length(sigma) == nrow(X))
                              k <- function(x, y) linear(x, y, sigma)
                              super$initialize(X, y, noise, k)
                            }
                          )
)

#' @export
GPR.polynomial <- R6::R6Class("GPR.polynomial", inherit = GPR,
                              public = list(
                                initialize = function(X, y, noise, sigma = fit(X,y,noise,"polynomial")$par[[1]], p = fit(X,y,noise,"polynomial")$par[[2]]){
                                  stopifnot(length(sigma) == 1, length(p) == 1)
                                  k <- function(x, y) polynomial(x, y, sigma, p)
                                  super$initialize(X, y, noise, k)
                                }
                              )
)

#' @export
GPR.sqrexp <-  R6::R6Class("GPR.sqrexp", inherit = GPR,
                           public = list(
                             initialize = function(X, y, noise, l = fit(X, y, noise, "sqrexp")$par){
                               stopifnot(length(l) == 1)
                               k <- function(x, y) sqrexp(x, y, l)
                               super$initialize(X, y, noise, k)
                             }
                             
                           )
)

#' @export
GPR.gammaexp <- R6::R6Class("GPR.gammaexp", inherit = GPR,
                          public = list(
                            initialize = function(X, y, noise, gamma = fit(X,y,noise,"gammaexp")$par[[1]], l = fit(X,y,noise,"gammaexp")$par[[2]]){
                              stopifnot(length(gamma) == 1, length(l) == 1)
                              k <- function(x, y) gammaexp(x, y, l, gamma)
                              super$initialize(X, y, noise, k)
                            }
                          )
)

#' @export
GPR.rationalquadratic <- R6::R6Class("GPR.rationalquadratic", inherit = GPR,
                            public = list(
                              initialize = function(X, y, noise, alpha = fit(X,y,noise,"rationalquadratic")$par[[1]], l = fit(X,y,noise,"rationalquadratic")$par[[2]]){
                                stopifnot(length(alpha) == 1, length(l) == 1)
                                k <- function(x, y) rationalquadratic(x, y, l, alpha)
                                super$initialize(X, y, noise, k)
                              }
                            )
)

# Funktion, um Kovarianzmatrix zu berechnen. Die verwendete Kovarianzfunktion muss bei Eingabe zweier Matrizen 
# gleicher Dimension die Werte bei Anwendung auf die jeweils i-ten Spalten für i = 1,...,ncol zurueckgeben.
covariance_matrix <- function(A, B, covariance_function) {
  outer(1:ncol(A), 1:ncol(B), function(i, j) covariance_function(A[, i, drop = F], B[, j, drop = F]))
}

#' @export
multivariate_normal <- function(n, mean, covariance, tol = 1e-6) {
  L <- tryCatch(error = function(cond) return(NULL), t(chol(covariance)))
  if (is.null(L)) {
    eig <- eigen(covariance, symmetric = TRUE)
    eigval <- eig$values
    stopifnot(all(eigval > - tol*abs(eigval[1])))
    L <- eig$vectors %*% diag(sqrt(pmax(eigval,0)))
  }
  drop(mean) + L %*% matrix(rnorm(n*length(mean), 0, 1), nrow = length(mean))
}

expand_range <- function(x) {
  r <- range(x)
  m <- mean(r)
  c(m - 1.2*(m - r[1]), m + 1.2*(m = r[2]))
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
rationalquadratic.matrix <- function(x, y, l, alpha) (1 + colSums((x - y)^2) / (2 * alpha * l^2))^(-alpha)
rationalquadratic.numeric <- function(x, y, l, alpha) (1 + sum((x - y)^2) / (2 * alpha * l^2))^(-alpha)