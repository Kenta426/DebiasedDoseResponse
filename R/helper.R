# Helper function for generating pseudo-outcomes
.pseudo.outcomes <- function(y, a, w, mu, g){
  n <- length(a)
  W.new <- data.frame(w[rep(1:n, n), ])
  w <- data.frame(w); colnames(w) <- colnames(W.new)
  A.new <- rep(a, each=n)
  muhat.mat <- matrix(mu(A.new, W.new), byrow=TRUE, nrow=n, ncol=n)
  muhat.obs <- diag(muhat.mat)
  mhat.obs <- rowMeans(muhat.mat)
  ghat.obs <- g(a, w)
  list(pseudo.outcome=(y - muhat.obs)/(ghat.obs) + mhat.obs,
       muhat.mat=muhat.mat,
       mhat.obs=mhat.obs)
}

# Generate kernel function
.kern = function(u, kernel="epa"){
  if (kernel=="epa") w <- 0.75*(1-u^2)*(abs(u)<=1)
  if (kernel=="uni") w <- 0.5*(abs(u)<=1)
  if (kernel=="tri") w <- (1-abs(u))*(abs(u)<=1)
  if (kernel=="gau") w <- dnorm(u)
  return(w)
}

# Matrix inverse via Cholesky decomposition
.qrXXinv = function(x, ...) {
  chol2inv(chol(crossprod(x)))
}

# Parse control argument for debiased_inference
.parse.debiased_inference <- function(...){
  option <- list(...); arg <- list()
  if (is.null(option$kernel.type)){
    kernel.type <- "epa"
  }
  else{
    kernel.type <- option$kernel.type
  }
  if (is.null(option$bandwidth.method)){
    bandwidth.method <- "DPI"
  }
  else{
    bandwidth.method <- option$bandwidth.method
  }
  if (is.null(option$alpha.pts)){
    alpha.pts <- 0.05
  }
  else{
    alpha.pts <- option$alpha.pts
  }
  if (is.null(option$unif)){
    unif <- TRUE
  }
  else{
    unif <- option$unif
  }
  if (is.null(option$alpha.unif)){
    alpha.unif <- 0.05
  }
  else{
    alpha.unif <- option$alpha.unif
  }
  if (is.null(option$bootstrap)){
    bootstrap <- 1e4
  }
  else{
    bootstrap <- option$bootstrap
  }
  arg$kernel.type <- kernel.type
  arg$bandwidth.method <- bandwidth.method
  arg$alpha.pts <- alpha.pts
  arg$unif <- unif
  arg$alpha.unif <- alpha.unif
  arg$bw.seq <- option$bw.seq
  arg$bootstrap <- bootstrap
  return(arg)
}
