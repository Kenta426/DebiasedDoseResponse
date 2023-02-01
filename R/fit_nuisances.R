.fit.regression <- function(Y, A, W,
                           SL.library = c("SL.mean", "SL.earth", "SL.gam",
                                          "SL.rpart", "SL.glm"),
                           cdf = TRUE){
  W <- data.frame(W)
  n <- length(A)
  if(cdf){
    X.ecdf <- ecdf(A); U <- X.ecdf(A);
  }else{
    U <- A
  }
  mu <- SuperLearner(Y=Y, X=cbind(U=U, W), SL.library = SL.library) # TODO: Automatically set to binomial for binary
  mu.hat <- function(x, w){
    w <- data.frame(w)
    names(w) <- names(w)
    if(cdf){
      return(c(predict(mu, newdata=cbind(U=X.ecdf(x), w), onlySL=TRUE)$pred))
    }else{
      return(c(predict(mu, newdata=cbind(U=x, w), onlySL=TRUE)$pred))
    }
  }
  return(mu.hat)
}

.cross.fit.regression <- function(Y, A, W,
                                  SL.library = c("SL.mean", "SL.earth", "SL.rpart",
                                                 "SL.gam", "SL.glm"),
                                  cdf = TRUE, folds = 5){
  n <- length(A)
  W <- data.frame(W)
  folds.id <- sample(rep(1:folds, length.out=n), size=n, replace=FALSE)
  muhats.fun.cv <- list() # store all functions
  X.ecdf.cv <- list() # store all functions

  for (v in 1:folds){
    train.inds <- which(folds.id != v)
    if(cdf){
      X.ecdf.cv[[v]] <- ecdf(A[train.inds]); U.v <- (X.ecdf.cv[[v]])(A[train.inds]);
    }else{
      U.v <- A[train.inds]
    }
    muhats.fun.cv[[v]] <- SuperLearner(Y=Y[train.inds],
                                       X=cbind(U=U.v, W[train.inds, ]),
                                       SL.library = SL.library)
  }

  mu.hat <- function(x, w){
    w <- data.frame(w)
    names(w) <- names(w)
    mu.preds <- rep(NA, length(x))
    if (length(x) != n){
      folds.id = rep(folds.id, n)
    }
    for (v in 1:folds){
      test.inds <- which(folds.id == v)
      if(cdf){
        mu.preds[test.inds] <- c(predict(muhats.fun.cv[[v]],
                                         newdata=cbind(U=(X.ecdf.cv[[v]])(x[test.inds]), w[test.inds, ]),
                                         onlySL=TRUE)$pred);
      }else{
        mu.preds[test.inds] <- c(predict(muhats.fun.cv[[v]],
                                         newdata=cbind(U=x[test.inds], w[test.inds, ]),
                                         onlySL=TRUE)$pred);
      }
    }
    return(mu.preds)
  }
  return(mu.hat)
}

.fit.density <- function(A, W, SL.library=c("SL.mean", "SL.earth", "SL.gam",
                                            "SL.rpart", "SL.glm"),
                         cdf=TRUE){
  W <- data.frame(W)
  names(W) <- names(W)
  n <- length(A)
  if(cdf){
    X.ecdf <- ecdf(A); U <- X.ecdf(A);
  }else{
    U <- A
  }
  g <- CD.SuperLearner(X=U, W=W, SL.library = SL.library, n.folds=10,
                       n.bins=2:floor(length(U)/50)) # TODO: Let people select folds and folds?
  g.hat <- function(x, w){
    w <- data.frame(w)
    names(w) <- names(w)
    if(cdf){
      return(c(predict.CD.SuperLearner(g, newdata=data.frame(cbind(X=X.ecdf(x), w)))$SL.density))
    }else{
      return(c(predict.CD.SuperLearner(g, newdata=data.frame(cbind(X=x, w)))$SL.density))
    }
  }
  return(g.hat)
}

.cross.fit.density <- function(A, W, SL.library=c("SL.mean", "SL.earth",
                                                  "SL.gam"),
                              cdf = TRUE, folds = 5){
  n <- length(A)
  folds.id <- sample(rep(1:folds, length.out=n), size=n, replace=FALSE)
  if(cdf){
    X.ecdf <- ecdf(A); U <- X.ecdf(A);
  }else{
    U <- A
  }

  ghats.fun.cv <- list() # store all functions
  X.ecdf.cv <- list() # store all functions

  for (v in 1:folds){
    train.inds <- which(folds.id != v)
    if(cdf){
      X.ecdf.cv[[v]] <- ecdf(A[train.inds]); U.v <- (X.ecdf.cv[[v]])(A[train.inds]);
    }else{
      U.v <- A[train.inds]
    }

    ghats.fun.cv[[v]] <- CD.SuperLearner(X=U.v, W=W[train.inds, ],
                                         SL.library = SL.library,
                                         n.folds=5, n.bins=2:6)
  }
  g.hat <- function(x, w){
    w <- data.frame(w)
    names(w) <- names(w)
    g.preds <- rep(NA, length(x))
    for (v in 1:folds){
      test.inds <- which(folds.id == v)
      if(cdf){
        g.preds[test.inds] <- c(predict.CD.SuperLearner(ghats.fun.cv[[v]],
                                                        newdata=cbind(X=(X.ecdf.cv[[v]])(x[test.inds]), w[test.inds, ]))$SL.density);
      }else{
        g.preds[test.inds] <- c(predict.CD.SuperLearner(ghats.fun.cv[[v]],
                                                        newdata=cbind(X=x[test.inds], w[test.inds, ]))$SL.density);
      }
    }
    return(g.preds)
  }
  return(g.hat)
}
