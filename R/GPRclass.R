#needs package R6

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
                       #K <- outer(1:ncol(X), 1:ncol(X), function(i,j) k(X[, factor(i)], X[, factor(j)]))
                       n <- ncol(X)
                       K <- matrix(0, nrow = n, ncol = n)
                       for (i in 1:n) {
                         for (j in 1:n) {
                           K[i, j] <-  k(X[, i], X[, j])
                         }
                       }
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
                     fit = function(){
                       K <- matrix(0, nrow = n, ncol = n)
                       for (i in 1:n) {
                         for (j in 1:n) {
                           K[i, j] <-  k(X[, i], X[, j])
                         }
                       }
                       private$.L <- t(chol(K + noise * diag(n)))
                       -0.5 * self$y %*%  - sum(log(diag(self$L))) - ncol(self$X) / 2 * log(2 * pi)
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
                              k <- function(x, y) c
                              super$initialize(X, y, k, noise)
                            }
                          )
)

GPR.linear <- R6::R6Class("GPR.linear", inherit = GPR,
                          public = list(
                            initialize = function(X, y, sigma, noise){
                              stopifnot(length(sigma) == nrow(X))
                              k <- function(x, y) sum(sigma * x * y)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

GPR.polynomial <- R6::R6Class("GPR.polynomial", inherit = GPR,
                              public = list(
                                initialize = function(X, y, sigma, p, noise){
                                  stopifnot(length(sigma) == 1, length(p) == 1)
                                  k <- function(x, y) (x %*% y + sigma)^p
                                  super$initialize(X, y, k, noise)
                                }
                              )
)

GPR.sqrexp <-  R6::R6Class("GPR.sqrexp", inherit = GPR,
                           public = list(
                             initialize = function(X, y, l, noise){
                               stopifnot(length(l) == 1)
                               k <- function(x, y,l=l) exp(-dist(rbind(x, y))^2/(2 * l^2))
                               super$initialize(X, y, k, noise)
                             }
                             
                           )
)

GPR.gammaexp <- R6::R6Class("GPR.gammaexp", inherit = GPR,
                          public = list(
                            initialize = function(X, y, gamma, l, noise){
                              stopifnot(length(gamma) == 1, length(l) == 1)
                              k <- function(x, y) exp(-(dist(rbind(x, y)) / l) ^ gamma)
                              super$initialize(X, y, k, noise)
                            }
                          )
)

GPR.rationalquadratic <- R6::R6Class("GPR.rationalquadratic", inherit = GPR,
                            public = list(
                              initialize = function(X, y, alpha, l, noise){
                                stopifnot(length(alpha) == 1, length(l) == 1)
                                k <- function(x, y) (1 + dist(rbind(x, y))^2 / (2 * alpha * l^2))^(-alpha)
                                super$initialize(X, y, k, noise)
                              }
                            )
)

cov_dict <- list(sqrexp = list(func = function(x, y,l) exp(-dist(rbind(x, y))^2/(2 * l^2)), deriv = function(x,y,l){
  r <- dist(rbind(x - y))
  r^2/l^3*exp(-r^2/(l^2*2))
}, start = c(1)
),
gammaexp = list(
  func = function(x, y, gamma, l) exp(-(dist(rbind(x, y)) / l) ^ gamma), deriv = function(x, y, gamma, l){
  r <- dist(rbind(x-y))
  c(-exp(-(r/l)^gamma) * (r/l)^gamma * log(r/l), exp(-(r/l)^gamma) * gamma * r^gamma / (l^(gamma+1)))
}, start = c(1,1)
),
constant = list(
  func = function(x, y, c) c, deriv = function(x, y, c) 0, start = c(1)
),
linear = list(
  func = function(x, y, sigma) sum(sigma * x * y), deriv = function(x, y, sigma) sigma, start = c(1)
),
polynomial = list(
  func = function(x, y, sigma, p) (x %*% y + sigma)^p, deriv = function(x, y, sigma, p){
    c(p * (x %*% y + sigma)^(p - 1), (x %*% y + sigma)^p * log((x %*% y + sigma)))
  }, start = c(1)
)
)

fit <-  function(X,y,noise,cov_names){
  param <- list()
  score <- c()
  for (cov in cov_names){
    usedcov <- cov_dict[[cov]]
    dens <- function(v){
      n <- ncol(X)
      K <- matrix(0, nrow = n, ncol = n)
      for (i in 1:n) {
        for (j in 1:n) {
          K[i, j] <-  do.call(usedcov$func, as.list(c(X[, i], X[, j],v)))
        }
      }
      
      L <- t(chol(K + noise * diag(n)))
      alpha <- solve(t(L), solve(L, y))
      - 0.5 * y %*% alpha - sum(log(diag(L))) - ncol(X) / 2 * log(2 * pi)
    }
    dens_deriv <- function(v){
      n <- ncol(X)
      K <- matrix(0, nrow = n, ncol = n)
      K_deriv <- array(0, c(n, n, length(usedcov$start)))
      for (i in 1:n) {
        for (j in 1:n) {
          K[i, j] <-  usedcov$func(X[, i], X[, j],...)
          K_deriv[i, j,] <- do.call(usedcov$deriv, as.list(c(X[, i], X[, j],v)))
        }
      }
      
      K_inv <- solve(K)
      alpha <- K_inv %*% y
      vapply(1:length(usedcov$start), function(i) 0.5 * sum(diag(alpha %*% t(alpha) - K_inv) %*% K_deriv[,,i]), numeric(1))
    }
    p <- optim(usedcov$start, dens, gr = dens_deriv, control = list(fnscale = -1))
    param <- append(param, p$par)
    print(param)
    score <- c(score, p$value)
  }
  return(list(par = param[[which.max(score)]],cov = cov_names[[which.max(score)]], score = score))
}


X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
noise <- 1
y <- c(0.1*X^3 + rnorm(length(X),0, 1))
kappa <- function(x,y) exp(-10*(x - y)^2)
Gaussian <- GPR$new(X, y, kappa, noise)
Gaussian$plot(seq(-5,5, by = 0.1))
z <- fit(X,y,noise,list("gammaexp"))
print(z)
Gaussian <- GPR$new(X, y, function(x, y) exp(-dist(rbind(x, y))^2/(2 * z$par^2)), noise)
Gaussian$plot(seq(-5,5, by = 0.1))
