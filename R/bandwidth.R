# Compute hat matrix for leave-one-out cross validation for a linear smoother
.hatmatrix <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
  }
  neval <- length(eval.pt)
  c.2 <- integrate(function(t){t^2*.kern(t, kernel.type)}, -Inf, Inf)$value
  e1 <- matrix(c(1, 0), ncol=2); e3 <- matrix(c(1, 0, 0, 0), ncol=4)
  w.h.zero   <- .kern(0, kernel.type)/h
  w.b.zero   <- .kern(0, kernel.type)/b
  hat.mat <- rep(0, neval)
  for (i in 1:neval){
    bw.min   <- sort(abs(x-eval.pt[i]))[10]
    h0     <- max(h, bw.min); b0     <- max(b, bw.min)
    k.h   <- .kern((x-eval.pt[i])/h0, kernel.type)/h0
    k.b   <- .kern((x-eval.pt[i])/b0, kernel.type)/b0
    ind.h <- k.h>0;  ind.b <- k.b>0
    N.h   <- sum(ind.h);  N.b <- sum(ind.b); ind   <- ind.b
    if (h>b) ind <- ind.h
    eY  <- y[ind]; eX  <- x[ind]
    K.h <- k.h[ind]; K.b <- k.b[ind]
    W.b <- matrix(NA, sum(ind), 4) # (1, u, u^2, u^3)
    for (j in 1:4)  W.b[,j] <- ((eX-eval.pt[i])/b0)^(j-1)
    W.h  <- matrix(NA, sum(ind), 2) # (1, u)
    for (j in 1:2)  W.h[,j] <- ((eX-eval.pt[i])/h0)^(j-1)
    invD.b  <- .qrXXinv((sqrt(K.b)*W.b))
    invD.h  <- .qrXXinv((sqrt(K.h)*W.h))
    hat.mat[i] <- (invD.h%*%t(e1*w.h.zero))[1,] -
      ((h/b)^2*c.2*invD.b%*%t(e3*w.b.zero))[3,]
  }
  approx(eval.pt, hat.mat, xout = x)$y
}


# Perform leave-one-out cross-validation based on the "hat matrix trick"
.robust.loocv <- function(x, y, h, b, eval.pt=NULL, kernel.type="epa"){
  ind <- order(x); x <- x[ind]; y <- y[ind]
  if (is.null(eval.pt)){
    eval.pt <- seq(quantile(x, 0.05), quantile(x, 0.95), length.out=30)
  }

  # Compute hat matrix --------------------------------------------------------
  hat.val <- .hatmatrix(x, y, h, b, eval.pt=eval.pt, kernel.type=kernel.type)

  # Compute debiased local linear estimator -----------------------------------
  est <- .lprobust(x, y, h, b, eval.pt=eval.pt, kernel.type=kernel.type)
  est.fn <- approx(eval.pt, est[,"theta.hat"], xout = x)$y

  # Compute the estimate of LOOCV risk ----------------------------------------
  mean(((y - est.fn) / (1 - hat.val))^2, na.rm=TRUE)
}
