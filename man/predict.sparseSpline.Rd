% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sparseSplines.R
\name{predict.sparseSpline}
\alias{predict.sparseSpline}
\title{Extract fitted values from a \code{sparseSpline} object}
\usage{
\method{predict}{sparseSpline}(object, newdata = NULL, which.lam = 10, ...)
}
\arguments{
\item{object}{an object of class \code{sparseSpline}, as produced by \code{\link{spam}}}

\item{newdata}{if NULL, returns predictions at original values, otherwise an N1 x p matrix with new X values at which to form predictions}

\item{which.lam}{an optional value for the index of the tuning parameter sequence to produce. If \code{NULL}, predicted values for all values of the tuning parameter.}

\item{...}{optional additional arguments for method compatability}
}
\value{
A matrix of dimension Nxd or an array of dimension N x nlambda x d if \code{which.lam=NULL}.
}
\description{
Extract fitted values from a \code{sparseSpline} object
}
\examples{
someData = generateSPASM(100,3,4)
SPAMfit = spam(someData$y, someData$x, returnall=TRUE)
filt = lgf1(someData$y, SPAMfit, someData$Tt,
  solve(someData$HHt), someData$GGt, rep(0,3), diag(1,3))
SPAMpreds = predict(SPAMfit, filt$xtt)
matplot(someData$y, ty='l', col=1, lty=1)
matlines(SPAMpreds, col=4, lty=3)
}
