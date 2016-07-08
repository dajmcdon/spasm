# spasm

The goal of `spasm` is to implement the first order Laplace Gaussian filter (as described in [Koyama et al](http://dx.doi.org/10.1198/jasa.2009.tm08326)) for estimating state space models whose observation equation has a mean function which is a spline in the unobserved state:
\[
\begin{align*}
y_t &= Z(x_t) + \epsilon_t\\
x_{t+1} &= T x_t + \eta_t,
\end{align*}
\]
where $Z(x)$ is a sparse spline estimated via the group lasso.

## Installation

You can install spasm from github with:

```R
# install.packages("devtools")
devtools::install_github("spasm/dajmcdon")
```

## Example

This is a basic example which shows you how to solve a common problem:

```R
someData = generateSPASM(100,3,4)
SPAMfit = spam(someData$y, someData$x)
filt = lgf1(someData$y, SPAMfit, someData$Tt,
  solve(someData$HHt), someData$GGt, rep(0,3), diag(1,3))
matplot(someData$x, ty='l', col=1, lty=1)
matlines(filt$xtt, col=2, lty=3)
```
