.kern = function(u, kernel="epa"){
  if (kernel=="epa") w <- 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w <- 0.5*(abs(u)<=1)
  if (kernel=="tri") w <- (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w <- dnorm(u)
  return(w)
}

.qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}

.lprobust <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)

  estimate.mat <- matrix(NA, neval, 4)
  colnames(estimate.mat) <- c("eval", "theta.hat", "mu.hat", "b.hat")
  c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value

  for (i in 1:neval){
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0     <- max(h, bw.min)
    b0     <- max(b, bw.min)

    k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0

    ind.h <- k.h>0;  ind.b <- k.b>0
    ind   <- ind.b
    if (h>b) ind <- ind.h
    eN  <- sum(ind); eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]

    W.b <- matrix(NA, eN, 4) # (1, u, u^2, u^3)
    for (j in 1:4)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2) # (1, u)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)

    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))

    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y
    beta.bias <- (h/b)^2*invD.b%*%crossprod(W.b*K.b, eY)*c.2
    tau <- beta.ll[1,1]
    tau.bias <- beta.bias[3,1]

    estimate.mat[i,] <- c(eval.pt[i], tau-tau.bias, tau, tau.bias)
  }
  estimate.mat
}

.locpoly <- function(x, y, h, eval.pt=NULL, kernel.type="epa", degree=1){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 2)
  colnames(estimate.mat) <- c("eval", "theta.hat")

  for (i in 1:neval){
    # prevents noninvertible D matrix
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0     <- max(h, bw.min)
    w.h   <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0

    ind <- w.h>0  # choose index for wider bandwith (more index). same if rho = 1
    eY  <- y[ind]; eX  <- x[ind]
    W.h <- w.h[ind] # kernelized X with h

    R.p <- matrix(NA, sum(ind), degree+1) # (1, u) for p = 1
    for (j in 1:(degree+1))  R.p[,j] <- ((eX-eval.pt[i])/h0)^(j-1)

    invG.p  <- .qrXXinv((sqrt(W.h)*R.p)) # (R.p^T W.h R.p)^-1
    beta.p  <- invG.p%*%crossprod(R.p*W.h, eY) #R.p^T W.h Y
    tau <- beta.p[1,1]
    estimate.mat[i,] <- c(eval.pt[i], tau)
  }
  estimate.mat
}

.loclinear <- function(x, y, h, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 2)
  colnames(estimate.mat) <- c("eval","theta.hat")

  # for each evaluation point
  for (i in 1:neval){
    # prevents noninvertible D matrix
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0     <- max(h, bw.min)
    w.h   <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0

    ind <- w.h>0  # choose index for wider bandwith (more index). same if rho = 1
    eY  <- y[ind]; eX  <- x[ind]
    W.h <- w.h[ind] # kernelized X with h

    R.p <- matrix(NA, sum(ind), 2) # (1, u) for p = 1
    for (j in 1:2)  R.p[,j] <- ((eX-eval.pt[i])/h0)^(j-1)

    invG.p  <- .qrXXinv((sqrt(W.h)*R.p)) # (R.p^T W.h R.p)^-1
    beta.p  <- invG.p%*%crossprod(R.p*W.h, eY) #R.p^T W.h Y
    tau <- beta.p[1,1]
    estimate.mat[i,] <- c(eval.pt[i], tau)
  }
  estimate.mat
}
