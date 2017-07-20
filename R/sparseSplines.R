
#' Estimates sparse (mulitivariate) additive model
#'
#' @param Y response, Nxd matrix
#' @param X predictors, Nxp matrix
#' @param df the number of spline basis functions to use, defaluts to 10
#' @param which.lam component of the tuning parameter vector for saving coefficients
#' @param returnall if TRUE, returns extra information, see below
#' @param ... optional arguments to \code{\link[SAM]{samQL}}
#'
#' @return a list with components
#' \describe{
#' \item{knots}{a list with components \code{interior}, a (df-1)xp matrix and \code{boundary}, a 2xd matrix}
#' \item{B.hat}{a (d)x(1+df*p) matrix containing the retained coefficients (intercept in the first column)}
#' \item{range.info}{a list with components \code{min} and \code{range} of size px1 each describing the predictors}
#' \item{spam.out}{if \code{returnall} is TRUE, a list of length d containing all output from \code{\link[SAM]{samQL}}. Useful for the \code{\link{predict}} or \code{\link{fitted}} methods}
#' \item{X}{the original \code{X} input}
#' }
#' @export
#'
#' @examples
#' someData = generateSPASM(100,3,4)
#' SPAMfit = spam(someData$y, someData$x)
spam <- function(Y, X, df = 10, which.lam = 10, returnall=FALSE,...){
  ## Requires SAM::samQL
  ## Inputs: Y - a N x d matrix,
  ##             X - a N x p,
  ##             df - the number of spline basis functions
  ##             which.lam - the column of the coef matrix to use (1:30)
  ##
  ## Output: knots - a list with components 'interior' as a (df-1) x p
  ##                          and boundary as a 2 x d
  ##              B.hat  -  a d x df*p matrix
  Y = as.matrix(Y)
  d = ncol(Y)
  p = ncol(X)
  B.hat = matrix(0, nrow=d, ncol=df*p) # for intercept
  spam.out = vector(mode='list',d)
  for(i in 1:d){
    ## samQL args: (X, y, p = 3, lambda = NULL,nlambda = NULL,
    ##              lambda.min.ratio=0.005,thol=1e-05,max.ite=1e+05)
    ## Notes: p in samQL refers to the number of basis functions (our 'df')
    spam.out[[i]] = SAM::samQL(X, Y[,i], p=df,...)
    B.hat[i,] = spam.out[[i]]$w[,which.lam]
    ## chosen the same for each component, need to choose this
  }
  # need the knots, but same for all d
  interior = spam.out[[1]]$nkots # typo in samQL
  boundary = spam.out[[1]]$Boundary.knots
  knots = list('interior'=interior,'boundary'=boundary)
  range.info = list('min'=spam.out[[1]]$X.min,'range'=spam.out[[1]]$X.ran)
  out = list('knots'=knots,'B.hat'=B.hat,'range.info'=range.info)
  if(returnall) {
    out$spamout = spam.out
    out$X = X
  }
  class(out) = 'sparseSpline'
  return(out)
}

#' Extract fitted values from a \code{sparseSpline} object
#'
#' @param object an object of class \code{sparseSpline}, as produced by \code{\link{spam}}
#' @param which.lam an optional value for the index of the tuning parameter sequence to produce. If \code{NULL}, fitted values for all values of the tuning parameter
#' @param ... optional additional arguments for method compatability
#'
#' @return A matrix of dimension Nxd or an array of dimension N x nlambda x d if \code{which.lam=NULL}.
#' @export
#'
#' @examples
#' someData = generateSPASM(100,3,4)
#' SPAMfit = spam(someData$y, someData$x, returnall=TRUE)
#' SPAMpreds = fitted(SPAMfit)
fitted.sparseSpline <- function(object, which.lam=10,...){
  predict(object, which.lam=which.lam,...)
}


#' Extract fitted values from a \code{sparseSpline} object
#'
#' @param object an object of class \code{sparseSpline}, as produced by \code{\link{spam}}
#' @param newdata if NULL, returns predictions at original values, otherwise an N1 x p matrix with new X values at which to form predictions
#' @param which.lam an optional value for the index of the tuning parameter sequence to produce. If \code{NULL}, predicted values for all values of the tuning parameter.
#' @param ... optional additional arguments for method compatability
#'
#' @return A matrix of dimension Nxd or an array of dimension N x nlambda x d if \code{which.lam=NULL}.
#' @export
#'
#' @examples
#' someData = generateSPASM(100,3,4)
#' SPAMfit = spam(someData$y, someData$x, returnall=TRUE)
#' filt = lgf1(someData$y, SPAMfit, someData$Tt,
#'   solve(someData$HHt), someData$GGt, rep(0,3), diag(1,3))
#' SPAMpreds = predict(SPAMfit, filt$xtt)
#' matplot(someData$y, ty='l', col=1, lty=1)
#' matlines(SPAMpreds, col=4, lty=3)
predict.sparseSpline <- function(object, newdata=NULL, which.lam=10,...){
  if(is.null(newdata)) newdata = object$X
  out = sapply(object$spamout, SAM::predict.samQL, newdata=newdata)
  out = simplify2array(out)
  if(!is.null(which.lam)) out = drop(out[,which.lam,])
  return(out)
}


