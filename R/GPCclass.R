#'  Predictions and Plots for Gaussian process classification
#'
#'  Implements a gaussian process for classification and gives tools 
#'  to predict and plot its values for given testpoints
#' 
#'
#' @usage \preformatted{GPC <- GPC$new(X, y, noise, cov_fun)
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
     initialize = function(X, y, epsilon, k){
       stopifnot(is.numeric(X), is.vector(y), is.numeric(y))
       stopifnot(is.numeric(epsilon), epsilon > 0, is.function(k))
       # Ist Input X ein Vektor, wird dieser als einzeilige Matrix behandelt.
       if (!is.matrix(X)) dim(X) <- c(1, length(X))
       stopifnot(length(y) == ncol(X))
       n <- length(y)
       K <- covariance_matrix(X, X, k)
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
     predict_class = function(X_star) {
       if (!is.matrix(X_star)) dim(X_star) <- c(1, length(X_star))
       P <- private$.sigmoid(self$f_hat)
       W <- P * (1 - P)
       K_star <- covariance_matrix(self$X, X_star, self$k)
       fs_bar <- t(K_star) %*% ((self$y + 1)/2 - P)
       v <- solve(self$L, (sqrt(W) * K_star))
       Vfs <- self$k(X_star, X_star) - colSums(v * v)
       sapply(1:length(Vfs), function(i) integrate(function(z) 
                private$.sigmoid(z) * dnorm(z, mean = fs_bar[i, 1], sd = Vfs[i]), -Inf, Inf)$value)
     },
     plot = function(limits = c(expand_range(self$X)[1], expand_range(self$X)[2]), length.out = 100L){
       if(nrow(self$X) == 1){
         testpoints <- seq(limits[1], limits[2] , length.out = length.out)
       } else if(nrow(self$X) == 2){
         s <- seq(limits[1], limits[2], length.out = length.out)
         testpoints <- matrix(c(rep(s, each = length.out), rep(s, length.out)), nrow = 2, byrow = TRUE)
       }
       if (nrow(self$X) > 2) {
         warning("Plot function not available for this dimension")
         return
       }
       if (is.vector(testpoints)) dim(testpoints) <- c(1, length(testpoints))
       predictions <- self$predict_class(testpoints)
       if(nrow(self$X) == 1){
         dat <- data.frame(x = t(testpoints), y = predictions)
         g <- ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y)) +
           ggplot2::theme_classic() +
           ggplot2::scale_y_continuous("output, p(y = 1| x)") +
           ggplot2::geom_tile(ggplot2::aes(x = x, y = y, fill = factor(as.integer(y > 0.5))), 
                              height = Inf, alpha = 0.5) +
           ggplot2::guides(fill = ggplot2::guide_legend(title = "Labels")) +
           ggplot2::scale_fill_manual(values = c("red", "blue")) +
           ggplot2::geom_line() +
           ggplot2::geom_point(data = data.frame(xpoints = c(self$X), ypoints = pmax(0,self$y)), 
                               mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = 4)) +
           ggplot2::scale_shape_identity()
       }
       else if(nrow(self$X) == 2){
         dat <- data.frame(x.1 = testpoints[1,], x.2 = testpoints[2,], y = 2*as.integer(predictions >= 0.5) - 1)
          g <- ggplot2::ggplot(dat, inherit.aes = F, ggplot2::aes(x = x.1, y = x.2, fill = factor(y))) +
           ggplot2::theme_classic() +
           ggplot2::scale_y_continuous(expression(x[2])) +
           ggplot2::scale_x_continuous(expression(x[1])) +
           ggplot2::geom_tile() + 
           ggplot2::scale_fill_manual(values = c("red", "blue")) +
           ggplot2::guides(fill = ggplot2::guide_legend(title = "Labels")) +
           ggplot2::geom_point(inherit.aes = F, data = data.frame(xpoints = c(self$X[1,]), ypoints = c(self$X[2,])), 
                               mapping = ggplot2::aes(x = xpoints, y = ypoints, shape = factor(self$y))) +
           ggplot2::scale_shape_manual(values = c(4, 2)) +
           ggplot2::guides(shape = ggplot2::guide_legend(title = "Testpoints")) +
           ggplot2::scale_color_manual(values = c("red", "blue"))
         }
       return(list(plot = g, pred = predictions))
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

"
#section for testing:
f <- function(x) (sum(abs(x)) > 2.5) - (!(sum(x) > 2.5))
limits <- matrix(c(-4, 4, -4, 4), nrow = 2, byrow = TRUE)
k <- function(x, y) sqrexp(x, y, 1)
#1-dim testpoints
testpoints <- seq(min(self$X),max(self$X), length.out = 150)
#2-dim testpoints
len <- 50
s <- seq(min(self$X),max(self$X),length.out = len)
testpoints <- matrix(c(rep(s,each = len), rep(s, len)), nrow = 2, byrow = TRUE)
"