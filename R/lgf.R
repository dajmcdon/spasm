#' Calculates the convolution of two Gaussian distributions
#'
#' @param x the mean of \eqn{x_{t|t}}
#' @param P the variance of \eqn{x_{t|t}}
#' @param Tt the state transition matrix
#' @param GGt the variance of the transition noise
#'
#' @details This function calculates
#' \deqn{p(x_{t+1} | y_{t}) = \int dx_{t} p(x_{t+1}| x_{t}) p(x_{t} | y_{1:t})}
#' assuming that both distributions on the right are Gaussian. This function also assumes that the mean of \eqn{x_{t+1|t}} is linear in \eqn{x_{t}}:
#' \deqn{p(x_{t+1} | x_{t}) = N(Tt x_{t|t}, GGt).}
#' The result is therefore:
#' \deqn{p(x_{t+1} | y_{t}) = N(Tt x_{t|t}, Tt P Tt' + GGt).}
#'
#' @return a list with components \code{x}, the mean of \eqn{x_{t+1|t}}, \code{P} the variance of \eqn{x_{t+1|t}} and \code{Pinv}, the inverse of \code{P}
#' @export
#'
#' @examples
#' posteriorPredictive(rnorm(3), diag(1,3), matrix(rnorm(9),3), diag(1,3))
posteriorPredictive <- function(x, P, Tt, GGt){
  ## Implements the "predict" step of the filter. Takes as input
  ## x_{t|t} and P_{t|t} along with system matrices Tt and GGt and
  ## returns x_{t+1|t} and P_{t+1|t}
  x.tplus1 = Tt %*% x
  P.tplus1= Tt %*% tcrossprod(P, Tt) + GGt
  Pinv = solve(P.tplus1)
  return(list(x=x.tplus1, P=P.tplus1, Pinv = Pinv))
}



#' Calculates negative \eqn{l(x_t)} up to an additive constant
#'
#' @param x p-vector of current states (to minimize over)
#' @param prev.state list containing mean and variance of \eqn{\hat{p}(x_{t-1} | y_{t-1})} (as produced by \code{\link{posteriorPredictive}}
#' @param y d-vector of observations at time t
#' @param out.spam list of fitted splines for the observation equation, as produced by \code{\link{spam}}
#' @param HHinv inverse of the observation error variance
#'
#' @return The negative of the (approximate) \eqn{l(x)}
#' @details This function calculates
#' \deqn{-l(x) = -log p(y_{t} | x_{t}) p(x_{t} | y_{1:t-1}).}
#' As the second distribution in the above product is approximated with a Gaussian, the result is also an approximation. Furthermore, as the first order LGF (\code{\link{lgf1}}) requires minimization over \eqn{-l(x)}, this function calculates it only to a constant of porportionality (independent of x)).
#' @export
negL <- function(x, prev.state, y, out.spam, HHinv){
  ## Calculates negative l(x_t) up to an additive constant
  ## See pseudoCode.pdf for notation.  IMPORTANT: RETURNS -ell(state)
  ## Input: x - p-vector of current states (to minimize over)
  ##           prev.state - list containing mean and variance of
  ##                       p.hat(state_{t-1}|y_{t-1}), (from
  ##                       posteriorPredictive)
  ##           y - d-vector of observables at time t
  ##           out.spam - list of fitted splines for the observation
  ##                              equation
  ##           HHinv - inverse of the observation error variance
  ## Output: -2 * l(x)
  ## For SPASM: these will be given by prev. Laplace Approx
  knots = out.spam$knots
  B.hat = out.spam$B.hat
  d = length(y)
  p = length(x)
  dfp = ncol(B.hat)
  df = dfp/p
  interior = t(t(knots$interior)*out.spam$range.info$range + out.spam$range.info$min)
  boundary = t(t(knots$boundary)*out.spam$range.info$range + out.spam$range.info$min)
  Phi.state = rep(0, dfp)
  for(j in 1:p){
    index.sel = (j - 1) * df + c(1:df)
    Phi.state[index.sel] = splines::ns(x[j], knots = interior[,j],
                                       Boundary.knots = boundary[,j])
  }
  Z.state = B.hat %*% Phi.state + out.spam$intercept
  obs = crossprod(Z.state, HHinv)
  observational.part = -2*obs %*% y + obs %*% Z.state
  pri = crossprod(x, prev.state$Pinv)
  prior.part = -2* pri %*% prev.state$x + pri %*% x
  return(observational.part + prior.part) #NOTE: This is 2 * negative log likelihood
}

#' Calculates the gradient of \eqn{-l(x)}
#'
#' @inheritParams negL
#'
#' @return The gradient of the (approximate) \eqn{-l'(x)}, (p-vector)
#' @details This function calculates
#' \deqn{-l'(x) = d/dx -log p(y_{t} | x_{t}) p(x_{t} | y_{1:t-1}).}
#' As the second distribution in the above product is approximated with a Gaussian, the result is also an approximation.
#' @export
gradNegL <- function(x, prev.state, y, out.spam, HHinv){
  knots = out.spam$knots
  B.hat = out.spam$B.hat
  d = length(y)
  p = length(x)
  dfp = ncol(B.hat)
  df = dfp/p
  interior = t(t(knots$interior)*out.spam$range.info$range + out.spam$range.info$min)
  boundary = t(t(knots$boundary)*out.spam$range.info$range + out.spam$range.info$min)
  Phi.state = rep(0,dfp)
  Phi.prime.state = matrix(0,nrow=dfp, ncol=p)
  for(k in 1:p){
    index.sel = (k - 1) * df + c(1:df)
    Phi.state[index.sel] = nsDeriv(x[k], derivs = 0, knots = interior[,k],
                                   Boundary.knots = boundary[,k])
    Phi.prime.state[index.sel,k] = nsDeriv(x[k], derivs = 1, knots = interior[,k],
                                           Boundary.knots = boundary[,k])
  }
  Z.prime.state = B.hat %*% Phi.prime.state
  Z.state = B.hat %*% Phi.state + out.spam$intercept
  observational.part = drop(crossprod(Z.prime.state, HHinv) %*% (y-Z.state))
  prior.part = drop(prev.state$Pinv %*% (prev.state$x - x))
  return(-2*(observational.part + prior.part)) ## we multiplied
    ## -2*l, so we do the same here.
}

#' Calculates the negative Hessian of \eqn{l(x)}
#'
#' @inheritParams negL
#'
#' @return The negative Hessian of the (approximate) \eqn{l'(x)}, (pxp matrix). Note that this is the negative of the Hessian of \eqn{l(x)}, not \eqn{-l(x)}.
#' @details This function calculates
#' \deqn{-l''(x) = - d^2/dx^2 log p(y_{t} | x_{t}) p(x_{t} | y_{1:t-1}).}
#' As the second distribution in the above product is approximated with a Gaussian, the result is also an approximation.
#' @export
negHessL <- function(x, prev.state, y, out.spam, HHinv){
  ## Calculates negative hessian of l(x_t) IMPORTANT: not hessian of
  ## -l(x_t), but - hessian of l(x_t)
  ## See pseudoCode.pdf for notation.
  ## Input: x - p-vector of current states (to minimize over)
  ##           prev.state - list containing mean and variance of
  ##                       p.hat(state_{t-1}|y_{t-1}), (from
  ##                       posteriorPredictive)
  ##           y - d-vector of observables at time t
  ##           out.spam - list of fitted splines for the observation
  ##                              equation
  ##           HHinv - inverse of the observation error variance
  ## Output: pxp-matrix -1*l''(x) evaluated at x
  knots = out.spam$knots
  B.hat = out.spam$B.hat
  d = length(y)
  p = length(x)
  dfp = ncol(B.hat)
  df = dfp/p
  interior = t(t(knots$interior)*out.spam$range.info$range + out.spam$range.info$min)
  boundary = t(t(knots$boundary)*out.spam$range.info$range + out.spam$range.info$min)
  Phi.state = rep(0,dfp)
  Phi.prime.state = matrix(0, nrow=dfp, ncol=p)
  Phi.prime.prime.state = matrix(0, nrow=dfp, ncol=p)
  for(j in 1:p){
    index.sel = (j - 1) * df + c(1:df)
    Phi.state[index.sel] = nsDeriv(x[j], derivs = 0, knots = interior[,j],
                                   Boundary.knots = boundary[,j])
    Phi.prime.state[index.sel,j] = nsDeriv(x[j], derivs = 1, knots = interior[,j],
                                           Boundary.knots = boundary[,j])
    Phi.prime.prime.state[index.sel,j] = nsDeriv(x[j], derivs = 2, knots = interior[,j],
                                                 Boundary.knots = boundary[,j])
  }
  Z.prime.prime.state = B.hat %*% Phi.prime.prime.state
  Z.prime.state = B.hat %*% Phi.prime.state
  Z.state = B.hat %*% Phi.state + out.spam$intercept
  ZppH = crossprod(Z.prime.prime.state, HHinv)
  term1 = diag(ZppH %*% y)
  term2 = diag(ZppH %*% Z.state)
  term3 = crossprod(Z.prime.state, HHinv) %*% Z.prime.state
  return(-2*(term1 - (term2 + term3) - prev.state$Pinv) )
}

#' First order LGF
#'
#' Given a sparse-spline representation of the observation equation, this function performs filtering to derive the conditional distributions of the filtered (\eqn{x_{t|t}}) and predicted (\eqn{x_{t|t-1}}) hidden states.
#'
#' @param Y the observed data, a Nxd matrix
#' @param out.spam the sparse additive spline model associated with the observation equation, as produced by \code{\link{spam}}
#' @param Tt the state transition matrix, pxp
#' @param HHinv the inverse of the observation noise variance matrix
#' @param GGt the transition noise variance matris
#' @param x0 prior mean of the state variables
#' @param P0 prior variance of the state variables
#'
#' @return A list with components \code{xt} and \code{xtt}, the predicted and filtered state means respectively and \code{Pt} and \code{Ptt}, the associated variances.
#' @export
#'
#' @examples
#' someData = generateSPASM(100,3,4)
#' SPAMfit = spam(someData$y, someData$x)
#' filt = lgf1(someData$y, SPAMfit, someData$Tt,
#'   solve(someData$HHt), someData$GGt, rep(0,3), diag(1,3))
#' matplot(someData$x, ty='l', col=1, lty=1)
#' matlines(filt$xtt, col=2, lty=3)
lgf1 <- function(Y, out.spam, Tt, HHinv, GGt, x0, P0){
  ##
  ##
  N = nrow(Y)
  p = ncol(Tt)
  xtt = matrix(0, nrow=N, ncol=p) # filtered state (x_{t|t})
  xt = matrix(0, nrow=N+1, ncol=p) # predicted state (x_{t+1|t})
  xt[1,] = x0
  Ptt = array(0, c(N, p, p))
  Pt = array(0, c(N+1, p, p))
  Pt[1,,] = P0
  pred.dist = list(x=x0, P = P0, Pinv = solve(P0))
  for(tt in 1:N) {
    ## cat(t,'\n')
    out.optim = optim(par=xt[tt,], fn=negL, gr=gradNegL,
                      prev.state = pred.dist, y = Y[tt,], out.spam = out.spam,
                      HHinv = HHinv,
                      method='BFGS')#, control=list(trace=6))#,abstol=1e-10))
    if(out.optim$convergence > 0) cat('convergence issue: ', out.optim$convergence)
    ## Filtering step
    xtt[tt,] = out.optim$par
    Ptt[tt,,] = solve(negHessL(xtt[tt,], pred.dist, Y[tt,], out.spam, HHinv))
    ## Prediction step
    pred.dist = posteriorPredictive(xtt[tt,], Ptt[tt,,],Tt, GGt)
    xt[tt+1,] = pred.dist$x
    Pt[tt+1,,] = pred.dist$P
  }
  return(list(xtt = xtt, xt = xt, Ptt = Ptt, Pt = Pt))
}

# The below has a number of issues:
# (1) k needs to be calculated fully, not proportionally, (same for l)
# (2) need to appropriately optimize each coordiante marginally


# neg.kj <- function(xj, x, j, prev.state, y, out.spam, HHinv, const=1e5){
#   negL = neg.ell(x, prev.state, y, out.spam, HHinv)
#   return(-log(xj+const)+negL)
# }
#
# deriv.neg.kj <- function(xj, x, j, prev.state, y, out.spam, HHinv, const=1e5){
#   ## Calculates derivative of -kj(x_t)
#   ## See pseudoCode.pdf for notation.  IMPORTANT: RETURNS derivative
#   ## of -ell(state)
#   ## Input: x - p-vector of current states (to minimize over)
#   ##           j - the coordinate to take derivatives of
#   ##           prev.state - list containing mean and variance of
#   ##                       p.hat(state_{t-1}|y_{t-1}), (from
#   ##                       posteriorPredictive)
#   ##           y - d-vector of observables at time t
#   ##           out.spam - list of fitted splines for the observation
#   ##                              equation
#   ##           HHinv - inverse of the observation error variance
#   ## Output: p-vector -1*l'(x)
#   knots = out.spam$knots
#   B.hat = out.spam$B.hat
#   d = length(y)
#   p = length(x)
#   dfp = ncol(B.hat)-1
#   df = dfp/p
#   interior = knots$interior[,j] * out.spam$range.info$range[j] +
#     out.spam$range.info$min[j]
#   boundary = knots$boundary[,j] * out.spam$range.info$range[j] +
#     out.spam$range.info$min[j]
#   index.sel = (j - 1) * df + c(1:df) + 1 # there's an intercept
#   Phi.state = nsDeriv(xj, derivs = 0, knots = interior, Boundary.knots = boundary)
#   Phi.prime.state = nsDeriv(xj, derivs = 1, knots = interior,
#                             Boundary.knots = boundary)
#   Z.prime.state = B.hat[,index.sel] %*% drop(Phi.prime.state)
#   Z.state = B.hat[,index.sel] %*% drop(Phi.state)
#   observational.part = drop(crossprod(Z.prime.state, HHinv) %*% (y-Z.state))
#   prior.part = drop(prev.state$Pinv[j,] %*% (prev.state$x - x))
#   negGradL = observational.part + prior.part
#   dkdx = 1/(xj+const)
#   return(-1*(dkdx+negGradL)) ## we multiplied
#       ## -2*l, so we do the same here.
# }
#
# det.neg.hess.k <- function(x, prev.state, y, out.spam, HHinv, const=1e5){
#   negHessL = neg.hess.ell(x, prev.state, y, out.spam, HHinv)
#   ddgdxdx = diag(-1/(x+const)^2)
#   ret = determinant(-1*ddgdxdx + negHessL)
#   return(ret$modulus)
# }
#
#
# lgf2 <- function(Y, out.spam, Tt, HHinv, GGt, x0, P0, const = 1e5){
#   N = nrow(Y)
#   p = ncol(Tt)
#   xtt = matrix(0, nrow=N, ncol=p) # filtered state (x_{t|t})
#   xt = matrix(0, nrow=N+1, ncol=p) # predicted state (x_{t+1|t})
#   xt[1,] = x0
#   Ptt = array(0, c(N, p, p))
#   Pt = array(0, c(N+1, p, p))
#   Pt[1,,] = P0
#   pred.dist = list(x=x0, P = P0, Pinv = solve(P0))
#   xbar = double(p)
#   for(tt in 1:N) {
#     ## cat(t,'\n')
#     for(j in 1:p){
#       out.optim.kj = optim(xt[tt,j], fn=neg.kj, gr=deriv.neg.kj, x=xt[tt,],
#                            j=j, prev.state = pred.dist, y = Y[tt,], out.spam = out.spam,
#                            HHinv = HHinv, const=const,
#                            method='BFGS', control=list(trace=6))#,abstol=1e-10))
#       if(out.optim.kj$convergence > 0) cat('convergence issue: ',out.optim.kj$convergence)
#       xbar[j] = out.optim.kj$par
#     }
#     out.optim.ell = optim(par=xt[tt,], fn=neg.ell, gr=grad.neg.ell,
#                           prev.state = pred.dist, y = Y[tt,], out.spam = out.spam,
#                           HHinv = HHinv,
#                           method='BFGS', control=list(trace=6))#,abstol=1e-10))
#     if(out.optim.ell$convergence > 0) cat('convergence issue: ', out.optim.ell$convergence)
#
#     ## Filtering step
#     xhat = out.optim.ell$par
#     Ptt[tt,,] = solve(neg.hess.ell(xhat, pred.dist, Y[tt,], out.spam, HHinv))
#     ## All on the log scale
#     neg.exp.ell = out.optim.ell$value # want the negative (in denominator)
#     exp.k = -neg.kj(xbar, xbar, NULL, pred.dist, Y[tt,], out.spam, HHinv, const) # this was negative
#     det.ell = determinant(Ptt[tt,,])$modulus # det of inverse
#     det.k = det.neg.hess.k(xbar, prev.state, Y[tt,], out.spam, HHinv, const=const)
#     xtt[,t] = exp(-det.k/2 + exp.k - det.ell/2 + neg.exp.ell) - const# exponentiate
#
#     ## Prediction step
#     pred.dist = posteriorPredictive(xtt[tt,], Ptt[tt,,],Tt, GGt)
#     xt[tt+1,] = pred.dist$x
#     Pt[tt+1,] = pred.dist$P
#   }
#   return(list(xtt = xtt, xt = xt, Ptt = Ptt, Pt = Pt))
# }
