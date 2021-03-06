% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GPCclass.R
\name{GPC}
\alias{GPC}
\title{Predictions and Plots for Gaussian process classification}
\usage{
\preformatted{ GPC <- GPC$new(X, y, cov_fun, epsilon)

 GPC$predict_class(X*)
 GPC$plot(limits, length.out)
 GPC$plot(X*) }
}
\description{
Implements a Gaussian process for classification and gives tools to predict
and plot its values for given test points
}
\section{Arguments}{

\describe{
 \item{\code{X}}{matrix of inputs}

 \item{\code{y}}{numeric vector of targets with values in \{-1,1\}}

 \item{\code{cov_fun}}{the chosen covariance function of the Gaussian process}

 \item{\code{epsilon}}{constant determining the threshold for the necessary
 optimization}

 \item{\code{X*}}{a numeric vector as the test input}

 \item{\code{limits}}{a numeric vector with lower and upper bounds for the plot}

 \item{\code{length.out}}{an integer indicating the number of points which are
 getting plotted}
 }
}

\section{Methods}{


 \code{$predict_class()} returns a numeric vector with the predicted
 posterior probabilities of class 1

 \code{$plot()} displays the results of the predict function for a number of
 points between limits or for all testpoints in a nice plot
}

\examples{
X <- matrix(seq(-1, 1, by = 0.1), nrow = 1)
y <- 2 * as.integer(X > 0) - 1
Gaussian_classifier <- GPC$new(X, y, 1e-5, function(x,y) exp(-3 * (x - y)^2))
Gaussian_classifier$plot()


}
\references{
Rasmussen, Carl E.; Williams, Christopher K. I. (2006).	Gaussian
 processes for machine learning
}
