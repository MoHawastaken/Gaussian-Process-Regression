#needs package R6

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
                       stopifnot(is.matrix(X), is.vector(y), is.numeric(y))
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
                         print(Pi)
                         W <- (1 - Pi)*Pi
                         L <- t(chol(diag(n) + (sqrt(W) %o% sqrt(W)) * K))
                         b <- W*f + (y + 1)/2 - Pi
                         intermediate <- solve(L, sqrt(W)*(K %*% b))
                         intermediate <- solve(t(L), intermediate)
                         a <- b - sqrt(W)*intermediate
                         f <- c(K %*% a)
                         print(a)
                         print(f)
                         objective <- -sum(a * f)/2 + sum(log(1/(1 + exp(-y * f))))
                         print(objective)
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
                       private$.X <- X
                       private$.y <- y
                       private$.k <- k
                     },
                     predict_class = function(Xs){
                       Pi <- 1/(1 + exp(-self$f_hat))
                       W <- Pi*(1 - Pi)
                       ks <- sapply(1:ncol(self$X), FUN = function(i) self$k(self$X[, i], Xs))
                       fs_bar <- sum(ks * ((self$y + 1)/2 - Pi))
                       v <- solve(self$L, (sqrt(W)*ks))
                       Vfs <- self$k(Xs, Xs) - sum(v * v)
                       hilfs_func <- function(z) 1/(1+exp(-z))*(1/sqrt(2*pi*Vfs))*exp(-(z-fs_bar)^2/(2*Vfs))
                       PIs_hat <- integrate(hilfs_func, -Inf, Inf)[1]
                       return(PIs_hat)
                     },
                     plot = function(testpoints){
                       dat <- data.frame(x = testpoints, 
                                         y = t(sapply(testpoints, function(x) self$predict(x)[1:2])))
                       ggplot2::ggplot(dat, ggplot2::aes(x = x, y = y.1)) +
                         ggplot2::theme_classic() +
                         ggplot2::scale_y_continuous("output, f(x)") +
                         ggplot2::geom_line() +
                         ggplot2::geom_ribbon(ggplot2::aes(ymin = y.1 - 2*sqrt(max(y.2,0)),
                                                           ymax = y.1 + 2*sqrt(max(y.2,0))), alpha = 0.2) +
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

X <- matrix(seq(-1,1,by = 0.2), nrow = 1)
y <- as.integer(X > 0)
kappa <- function(x,y) exp(-(x - y)^2)
gaussian_classifier <- GPC$new(X, y, kappa, 0.1)