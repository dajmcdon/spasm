nsDeriv <- function(x, derivs = 0, df = NULL, knots = NULL, intercept = FALSE,
                    Boundary.knots = range(x)){
  ## This function is the same as splines::ns in base R
  ## The only change is to allow an argument which evaluates the derivatives.
  nx <- names(x)
  x <- as.vector(x)
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  if (!missing(Boundary.knots)) {
    Boundary.knots <- sort(Boundary.knots)
    outside <- (ol <- x < Boundary.knots[1L]) | (or <- x >
                                                   Boundary.knots[2L])
  }
  else outside <- FALSE
  ## my addition
  if(derivs > 3) stop('Too many derivatives requested, max is 3 for natural cubic spline')
  derivs = rep(derivs, length(x))
  ## end my addition
  if (!missing(df) && missing(knots)) {
    nIknots <- df - 1 - intercept
    if (nIknots < 0) {
      nIknots <- 0
      warning("'df' was too small; have used ", 1 + intercept)
    }
    knots <- if (nIknots > 0) {
      knots <- seq.int(0, 1, length.out = nIknots + 2L)[-c(1L,
                                                           nIknots + 2L)]
      stats::quantile(x[!outside], knots)
    }
  }
  else nIknots <- length(knots)
  ## cat('nIknots',nIknots)
  Aknots <- sort(c(rep(Boundary.knots, 4L), knots))
  ##  cat('Aknots')
  ##  cat(Aknots)
  if (any(outside)){
    ## warning('outside')
    basis <- array(0, c(length(x), nIknots + 4L))
    if (any(ol)) {
      k.pivot <- Boundary.knots[1L]
      xl <- cbind(1, x[ol] - k.pivot)
      tt <- splines::splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,1))
      basis[ol, ] <- xl %*% tt
    }
    if (any(or)) {
      k.pivot <- Boundary.knots[2L]
      xr <- cbind(1, x[or] - k.pivot)
      tt <- splines::splineDesign(Aknots, rep(k.pivot, 2L), 4, c(0,1))
      basis[or, ] <- xr %*% tt
    }
    if (any(inside <- !outside)) basis[inside, ] <- splines::splineDesign(Aknots, x[inside], 4)
  }
  else basis <- splines::splineDesign(Aknots, x, ord=4, derivs=derivs) # this line gets the derivatives
  const <- splines::splineDesign(Aknots, Boundary.knots, 4, c(2, 2))
  if (!intercept) {
    const <- const[, -1, drop = FALSE]
    basis <- basis[, -1, drop = FALSE]
  }
  qr.const <- qr(t(const))
  basis <- as.matrix((t(qr.qty(qr.const, t(basis))))[, -(1L:2L),drop = FALSE])
  n.col <- ncol(basis)
  if (nas) {
    nmat <- matrix(NA, length(nax), n.col)
    nmat[!nax, ] <- basis
    basis <- nmat
  }
  dimnames(basis) <- list(nx, 1L:n.col)
  a <- list(degree = 3, knots = if (is.null(knots)) numeric() else knots,
            Boundary.knots = Boundary.knots, intercept = intercept)
  attributes(basis) <- c(attributes(basis), a)
  class(basis) <- c("ns", "basis", "matrix")
  basis
}
