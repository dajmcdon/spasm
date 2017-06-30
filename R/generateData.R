
#' Generate data from an additive state space model
#'
#' This function generates data and returns components so as to easily test the Laplace-Gaussian filter and the sparse spline additive models. Data is generated as
#' \deqn{y_{t} = Z sin(z * x_{t} / 2\pi) + \sigma_{y} * N(0, I)}
#' \deqn{x_{t+1} = Tt x_t + \sigma_x * N(0,I)}
#' where \eqn{Z} is a dxp matrix of standard uniform random variables and \eqn{z=[1 ... p]}.
#'
#'
#' @param N the number of observed time points
#' @param p the dimension of the hidden state vector
#' @param d the dimension of the observation vector
#' @param sig.x scalar which multiplies the hidden standard normal noise
#' @param sig.y scalar which multiplies the standard normal noise added to the observations
#'
#' @return a list with components
#' \describe{
#' \item{y}{the observations, nxd}
#' \item{x}{the hidden states, nxp}
#' \item{Tt}{the state transition matrix, pxp}
#' \item{HHt}{the observation noise variance, dxd}
#' \item{GGt}{the hidden noise variance, pxp}
#' }
#' @export
#'
#' @examples
#' generateSPASM(100,3,4)
generateSPASM <- function(N, p, d, sig.x=.1, sig.y = .1){
 Tt = diag(.9, p)
 x = matrix(0, p, N)
 Zt = matrix(runif(p*d), d, p)
 for(i in 2:N){
   x[,i] = Tt %*% x[,i-1] + sig.x * rnorm(p)
 }
 y = Zt %*% sin((1:p)/(2*pi) * x) + sig.y * matrix(rnorm(d*N), d, N)
 return(list(y=t(y), x=t(x), Tt=Tt, HHt = diag(sig.y, d), GGt = diag(sig.x, p), Zt = Zt))
}
