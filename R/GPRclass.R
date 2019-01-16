#needs package R6

GPR <- R6::R6Class("GPR",
                   private = list(
                     .X = NA,
                     .k = NA,
                     .y = NA,
                     .L = NA,
                     .alpha = NA,
                     .noise = NA
                   ),
                   public = list(
                     initialize = function(X, y, k, noise){
                       stopifnot(is.matrix(X), is.vector(y), is.numeric(y))
                       stopifnot(is.numeric(noise), length(noise) == 1, is.function(k))
                       private$.X <-  X
                       private$.y <- y
                       private$.k <-  k
                       private$.noise <- noise
                       #K <- outer(1:ncol(X), 1:ncol(X), function(i,j) k(X[, factor(i)], X[, factor(j)]))
                       n = ncol(X)
                       K <- matrix(0, nrow = n, ncol = n)
                       for (i in 1:n) {
                         for (j in 1:n) {
                           K[i, j] = k(X[, i], X[, j])
                         }
                       }
                       private$.L <- chol(K + noise * diag(n))
                       private$.alpha <- solve(t(self$L), solve(self$L, y))
                     },
                      predict = function(Xs){
                       
                       #Calculate k* (ks)
                       ks <- sapply(1:ncol(self$X), FUN = function(i) self$k(self$X[, i], Xs))
                       
                       #calculate all other variables directly
                       fs <- ks %*% self$alpha
                       v <- solve(self$L, ks)
                       Vfs <- self$k(Xs, Xs) - v %*% v
                       logp <- -0.5 * self$y %*% self$alpha - sum(log(diag(self$L))) - ncol(self$X) / 2 * log(2 * pi)
                       return(c(fs, Vfs, logp))
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
                     }
                   )
  
)

X <- matrix((1:200)/10, nrow = 1)
y <- c(10*X^2)
kappa <- function(x,y) exp(-(x-y)^2)
noise <- 0.1
Gaussian <- GPR$new(X, y, kappa, noise)
plot(Vectorize(function(x) Gaussian$predict(x)[1]), 0, 20)

GPR.constant <- R6::R6Class("GPR.constant",
                          inherit = GPR,
                          public = list(
                            initialize = function(X, y, c, noise){
                              k <- function(x, y) c
                              super$initialize(X, y, k, noise)
                            }
                          )
)                          