#' Optimization of hyperparameters
#' 
#' Applies optimization methods to find optimal parameters of the covariance function for a given set of datapoints
#' 
#' @section Usage: 
#' \preformatted{fit(X, y, noise, cov_names)}
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
#'   \code{cov_names} a list of names of covariance functions
#' @name fit
#' @references Rasmussen, Carl E.; Williams, Christopher K. I. (2006).	Gaussian processes for machine learning
NULL

#Save derivatives of covariance functions
cov_dict <- list(
  sqrexp = list(func = sqrexp, display = "Squared Exponential",
                deriv = function(x, y, l){
                  r <- sqrt(sum((x - y)^2))
                  r^2/l^3*exp(-r^2/(l^2*2))
                }, start = c(1)
  ),
  gammaexp = list(func = gammaexp, display = "Gamma Exponential",
                  deriv = function(x, y, gamma, l){
                    r <- sqrt(sum((x - y)^2))
                    c(-exp(-(r/l)^gamma) * (r/l)^gamma * log(r/l), exp(-(r/l)^gamma) * gamma * r^gamma / (l^(gamma + 1)))
                  }, start = c(1, 1)
  ),
  constant = list(func = constant, display = "Constant", deriv = function(x, y, c) 0, start = c(1)
  ),
  linear = list(func = linear, display = "Linear", deriv = function(x, y, sigma) sigma, start = c(1)
  ),
  polynomial = list(func = polynomial, display = "Polynomial",
                    deriv = function(x, y, sigma, p){
                      c(p * (x %*% y + sigma)^(p - 1), (x %*% y + sigma)^p * log((x %*% y + sigma)))
                    }, start = c(1, 2)
  ),
  rationalquadratic = list(func = rationalquadratic, display = "Rational Quadratic",
                           deriv = function(x, y, alpha, l){
                              if (is.vector(x)) r <- sum((x - y)^2)
                              else r <- colSums((x - y)^2)
                              c(((r/(2*l^2*alpha) + 1)^(-alpha) * (r - (2*l^2*alpha + r) * log(r/(2 * l^2*alpha) 
                                                                                               + 1)))/(2*l^2*alpha + r),
                                (r*(r/(2*l^2*alpha) + 1)^(-alpha - 1))/(l^3))
                           }, start = c(1,1)
  )
)
cov_df <- data.frame(list(name = "sqrexp", display = "Squared Exponential", start = I(list(1)), 
                          func = I(list(cov_dict[["sqrexp"]]$func)), deriv = I(list(cov_dict[["sqrexp"]]$deriv)))
                     , stringsAsFactors = FALSE )

for (name in names(cov_dict)[-1]){
  d1 <- list(name = name, display = cov_dict[[name]]$display, 
                        start = I(list(cov_dict[[name]]$start)),
             func =  I(list(cov_dict[[name]]$func)), 
             deriv = I(list(cov_dict[[name]]$deriv)))
  cov_df[nrow(cov_df) + 1,] <-  I(d1)
}
row.names(cov_df) <- cov_df$name

#' @export
fit <-  function(X, y, noise, cov_names){
  param <- list()
  score <- c()
  for (cov in cov_names){
    #if(class(try(solve(K + noise * diag(n)),silent=T)) != "matrix"){
    #  warning("K(X,X) + noise * I is not invertible for ", cov,". The algorithm is not defined for this case.")
    #  param <- append(param, -Inf)
    #  score <- c(score, -Inf)
    #  next
    #}
    usedcov <- cov_df[cov,]
    nparam <- length(usedcov$start[[1]])
    l <- list() #parameters for optim()
    best_par <- list(usedcov$start[[1]])
    prev_par <- c()
    dens <- function(v){
      K <- covariance_matrix(X, X, function(x, y) do.call(usedcov$func[[1]], append(list(x, y), v)))

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
            K[i, j] <-  do.call(usedcov$func[[1]], append(list(X[, i], X[, j]), v))
            K_deriv[i, j,] <- do.call(usedcov$deriv[[1]], append(list(X[, i], X[, j]), v))
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
    if (cov == "polynomial"){
      p_s <- list()
      p_sc <- list()
      for (i in 1:10){
        q <- optim(usedcov$start[[1]][1], function(sig) dens(c(sig,i)), method = "Brent", 
              lower = 0, upper = 5, control = list(fnscale = -1))
        p_s <- append(p_s, list(q))
        p_sc <- append(p_sc, q$value)
        
      }
      p <- list(par = c(p_s[[which.max(unlist(p_sc))]]$par, which.max(unlist(p_sc))), value = max(unlist(p_sc)))
    }
    else{
    l <- append(l, list(control = list(fnscale = -1)))
    p <- do.call(optim, append(list(usedcov$start[[1]], dens), l))
    }
    print(p)
    param <- append(param, list(p$par))
    score <- c(score, p$value)
    
  }
  return(list(par = param[[which.max(score)]], cov = cov_names[[which.max(score)]], score = score))
}

#section for testing:

X <- matrix(seq(-5,5,by = 0.2), nrow = 1)
y <- c(0.1*X^3 + rnorm(length(X), 0, 1))

z <- fit(X, y, noise = 1, cov_names = list("polynomial", "sqrexp", "gammaexp"))

print(z)
Gaussian <- GPR$new(X, y, function(x,y) do.call(cov_df[z$cov, ]$func[[1]], append(list(x, y), z$par)), noise = 1)
Gaussian$plot(seq(-5, 5, by = 0.1))
