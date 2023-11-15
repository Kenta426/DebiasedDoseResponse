# Debiased local linear regression
.lprobust <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
  }
  neval <- length(eval.pt)
  estimate.mat <- matrix(NA, neval, 4)
  colnames(estimate.mat) <- c("eval", "theta.hat", "mu.hat", "b.hat")
  c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
  for (i in 1:neval){
    bw.min <- sort(abs(x-eval.pt[i]))[10]
    h0 <- max(h, bw.min); b0 <- max(b, bw.min)
    k.h <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0; ind <- ind.b
    if (h>b) ind <- ind.h
    eN <- sum(ind); eY <- y[ind]; eX <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, eN, 4)
    for (j in 1:4)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y
    beta.bias <- (h/b)^2*invD.b%*%crossprod(W.b*K.b, eY)*c.2
    estimate.mat[i,] <- c(eval.pt[i], beta.ll[1,1]-beta.bias[3,1],
                          beta.ll[1,1], beta.bias[3,1])
  }
  estimate.mat
}
.lprobust.quad <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]; n<-length(y)
  if (is.null(eval.pt)){
    x.max <- max(x); x.min <- min(x)
    eval.pt <- seq(x.min, x.max, length.out=30)
  }
  neval <- length(eval.pt)
  Estimate <- matrix(NA, neval, 4)
  colnames(Estimate) <- c("eval", "theta.hat", "mu.hat", "b.hat")
  for (i in 1:neval){
    # this part of code is simply avoiding a singular matrix
    bw.min   <- sort(abs(x-eval.pt[i]))[4]
    h0 <- max(h, bw.min); b0 <- max(b, bw.min)
    k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0; ind <- ind.b
    if (h>b) ind <- ind.h
    eN  <- sum(ind); eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, eN, 3) # (1, u, u^2)
    for (j in 1:3)  W.b[, j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h <- matrix(NA, eN, 2) # (1, u)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    # estimate fixed-h c.2
    u <- (eX-eval.pt[i])/h0
    c2.est <- (invD.h %*% crossprod(W.h*K.h,u^2))[1]
    beta.ll  <- invD.h%*%crossprod(W.h*K.h, eY) #R.p^T W.h Y
    beta.bias <- (h0/b0)^2*invD.b%*%crossprod(W.b*K.b, eY)*c2.est
    tau <- beta.ll[1,1] # local linear term
    tau.bias <- beta.bias[3,1] # estimated bias
    Estimate[i,] <- c(eval.pt[i], tau-tau.bias, tau, tau.bias)
  }
  Estimate
}


# # Local polynomial regression
# .locpoly <- function(x, y, h, eval.pt=NULL, kernel.type="epa", degree=1){
#   ind <- order(x); x <- x[ind]; y <- y[ind]
#   if (is.null(eval.pt)){
#     eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
#   }
#   neval <- length(eval.pt)
#   estimate.mat <- matrix(NA, neval, 2)
#   colnames(estimate.mat) <- c("eval", "theta.hat")
#   # Compute for each evaluation point
#   for (i in 1:neval){
#     bw.min   <- sort(abs(x-eval.pt[i]))[10]
#     h0     <- max(h, bw.min)
#     k.h   <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0
#     ind <- k.h>0
#     eY  <- y[ind]; eX  <- x[ind]; K.h <- k.h[ind]
#     W.h <- matrix(NA, sum(ind), degree+1)
#     for (j in 1:(degree+1))  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
#     invD.h <- .qrXXinv((sqrt(K.h)*W.h))
#     beta <- invD.h%*%crossprod(W.h*K.h, eY)
#     estimate.mat[i,] <- c(eval.pt[i], beta[1,1])
#   }
#   estimate.mat
# }
#
# # Local linear regression
# .loclinear <- function(x, y, h, eval.pt=NULL, kernel.type="epa"){
#   ind <- order(x); x <- x[ind]; y <- y[ind]
#   if (is.null(eval.pt)){
#     eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95),  length.out=30)
#   }
#   neval <- length(eval.pt)
#   estimate.mat <- matrix(NA, neval, 2)
#   colnames(estimate.mat) <- c("eval","theta.hat")
#   for (i in 1:neval){
#     bw.min   <- sort(abs(x-eval.pt[i]))[10]
#     h0 <- max(h, bw.min)
#     k.h <- .kern((x-eval.pt[i])/h0, kernel=kernel.type)/h0
#     ind <- k.h>0
#     eY <- y[ind]; eX <- x[ind]; K.h <- k.h[ind]
#     w.h <- matrix(NA, sum(ind), 2)
#     for (j in 1:2)  w.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
#     invD.h  <- .qrXXinv((sqrt(K.h)*w.h))
#     beta.p  <- invD.h%*%crossprod(w.h*K.h, eY)
#     estimate.mat[i,] <- c(eval.pt[i], beta.p[1,1])
#   }
#   estimate.mat
# }
