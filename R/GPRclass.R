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
                       private$.L <- chol(K + sigma * diag(n))
                       private$.alpha <- solve(t(L), solve(L, y))
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

GPR.constant <- R6::R6Class("GPR.constant",
                          inherit = GPR,
                          public = list(
                            initialize = function(X, y, c, noise){
                              k <- function(x, y) c
                              super$initialize(X, y, k, noise)
                            }
                          )
)                          