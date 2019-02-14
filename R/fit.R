#' Optimization of hyperparameters
#' 
#' Applies optimization methods to find optimal parameters of the covariance function for a given set of datapoints
#' 
#' @section Usage: 
#' \preformatted{fit(X, y, noise, cov_list)}
#'
#'
#' @section Arguments:
#' 
#'   \code{X} matrix of inputs
#'
#'   \code{y} numeric vector of targets
#' 
#'   \code{noise} the inflicted noise of the observations
#' 
#'   \code{cov_list} a list of names of covariance functions
#' @name fit
#' @references Rasmussen, Carl E. W., Christopher K. I. (2006).	Gaussian processes for machine learning

#Save derivatives of covariance functions for the optimization of their hyperparameters
cov_dict <- list(
  sqrexp = list(func = sqrexp, 
                deriv = function(x, y, l){
                  r <- sqrt(sum((x - y)^2))
                  r^2/l^3*exp(-r^2/(l^2*2))
                }, start = c(1)
  ),
  gammaexp = list(func = gammaexp, 
                  deriv = function(x, y, gamma, l){
                    r <- sqrt(sum((x - y)^2))
                    c(-exp(-(r/l)^gamma) * (r/l)^gamma * log(r/l), exp(-(r/l)^gamma) * gamma * r^gamma / (l^(gamma + 1)))
                  }, start = c(1, 1)
  ),
  constant = list(func = constant, deriv = function(x, y, c) 0, start = c(1)
  ),
  linear = list(func = linear, deriv = function(x, y, sigma) sigma, start = c(1)
  ),
  polynomial = list(func = polynomial, 
                    deriv = function(x, y, sigma, p){
                      c(p * (x %*% y + sigma)^(p - 1), (x %*% y + sigma)^p * log((x %*% y + sigma)))
                    }, start = c(1, 2)
  ),
  rationalquadratic = list(func = rationalquadratic,
                           deriv = function(x, y, alpha, l){
                              r <- sum((x - y)^2)
                              c(((r/(2*l^2*alpha) + 1)^(-alpha) * (r - (2*l^2*alpha + r) * log(r/(2 * l^2*alpha) 
                                                                                               + 1)))/(2*l^2*alpha + r),
                                (r*(r/(2*l^2*alpha) + 1)^(-alpha - 1))/(l^3))
                           }, start = c(1,1)
  )
)
#' @export
fit <-  function(X, y, noise, cov_names){
  param <- list()
  score <- c()
  for (cov in cov_names){
    usedcov <- cov_dict[[cov]]
    nparam <- length(usedcov$start)
    l <- list() #parameters for optimization
    dens <- function(v){
      K <- covariance_matrix(X, X, function(x,y) do.call(usedcov$func, append(list(x, y), v)))
      
      L <- t(chol(K + noise * diag(ncol(X))))
      alpha <- solve(t(L), solve(L, y))
      - 0.5 * y %*% alpha - sum(log(diag(L))) - ncol(X) / 2 * log(2 * pi)
    }
    if (cov %in% c("sqrexp", "gammaexp","rationalquadratic")){
      dens_deriv <- function(v){
        n <- ncol(X)
        K <- matrix(0, nrow = n, ncol = n)
        K_deriv <- array(0, c(n, n, nparam))
        for (i in 1:n) {
          for (j in 1:n) {
            K[i, j] <-  do.call(usedcov$func, as.list(c(X[, i], X[, j], v)))
            K_deriv[i, j,] <- do.call(usedcov$deriv, as.list(c(X[, i], X[, j], v)))
          }
        }
        
        K_inv <- solve(K)
        alpha <- K_inv %*% y
        vapply(1:nparam, function(i) 0.5 * sum(diag(alpha %*% t(alpha) - K_inv) %*% K_deriv[,,i]), numeric(1))
      }
      l <- append(l, list(gr = dens_deriv))
    }
    #Switch optim method if parameter is one dimensional
    if (nparam == 1) l <- append(l, list(method = "Brent", lower = 0, upper = 10))
    else l <- append(l, list(method = "BFGS"))
    l <- append(l, list(control = list(fnscale = -1)))
    p <- do.call(optim, append(list(usedcov$start, dens), l))
    param <- append(param, list(p$par))
    score <- c(score, p$value)
  }
  return(list(par = param[[which.max(score)]], cov = cov_names[[which.max(score)]], score = score))
}


X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
noise <- 1
y <- c(0.1*X^3 + rnorm(length(X),0, 1))

z <- fit(X,y,noise,list("rationalquadratic"))

print(z)
Gaussian <- GPR$new(X, y, function(x,y) do.call(cov_dict[[z$cov]]$func, append(list(x,y),z$par)), noise)
Gaussian$plot(seq(-5,5, by = 0.1))
