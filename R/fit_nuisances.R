# File: influence_function.R
# Author: Kenta Takatsu
# Description:
#   Two functions in this file corresponds to the influence functions
#   for debiased local linear estimators and local linear estimators.
# References:
#   Takatsu K., and Westling T., (2023).
#   "Debiased inference for a covariate-adjusted regression function"
#   Kennedy EH., Ma Z., McHugh MD, Small DS (2017).
#   "Nonparametric methods for doubly robust estimation of
#   continuous treatment effects"
# Diaz Munoz, Ivan, and Mark J. van der Laan. (2011)
#   "Super learner based conditional density
#   estimation with application to marginal structural models."

#' Title: .fit.regression
#' Run outcome regression based on SuperLearner package and returns new
#' function that generates values for new points
#'
#' @param Y \code{n x 1} numeric vector of observed outcome values.
#' @param A \code{n x 1} numeric vector of observed exposure values.
#' @param W \code{n x p} matrix of potential confounders.
#' @param SL.library list of base learners for SuperLearner package
#' @param cdf TRUE if it transforms the exposure A to\code{[0,1]} by empirical CDF
#'  the default is TRUE.
#' @return
#' @export
#' @examples
#' # Sample data
#' n <- 200; cols <- 3
#' W <- matrix(runif(n*cols), ncol=cols) # 200 * 3 matrix of covariates
#' A <- rnorm(n, mean=W%*%rnorm(cols))   # 200 * 1 vector of treatment variable
#' Y <- rnorm(n, mean=sin(A)+W%*%rnorm(cols))  # 200 * 1 vector of response
#' mu.hat <- .fit.regression(Y, A, W)
#' mu.hat(A, W) # estimated values
.fit.regression <- function(Y, A, W,
                            SL.library=c("SL.earth","SL.gam",
                                         "SL.glm","SL.glmnet",
                                         "SL.glm.interaction",
                                         "SL.mean"),
                           cdf = TRUE){
  W <- data.frame(W)
  n <- length(A)
  if(cdf){
    X.ecdf <- ecdf(A); U <- X.ecdf(A);
  }else{
    U <- A
  }
  mu <- SuperLearner(Y=Y, X=cbind(U=U, W), SL.library = SL.library)
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

#' Title: .fit.semiparametric.density
#' Estimate standardized propensity score using semiparametric model.
#' The implementation is based on npcausal package by Edward Kennedy.
#'
#' @param A \code{n x 1} numeric vector of observed exposure values.
#' @param W \code{n x p} matrix of potential confounders.
#' @param SL.library list of base learners for SuperLearner package
#' @return
#' @export
#' @examples
#' # Sample data
#' n <- 200; cols <- 3
#' W <- matrix(runif(n*cols), ncol=cols) # 200 * 3 matrix of covariates
#' A <- rnorm(n, mean=W%*%rnorm(cols))   # 200 * 1 vector of treatment variable
#' Y <- rnorm(n, mean=sin(A)+W%*%rnorm(cols))  # 200 * 1 vector of response
#' mu.hat <- .fit.semiparametric.density(A, W)
#' mu.hat(A, W) # estimated values
.fit.semiparametric.density <- function(A, W,
                                        SL.library=c("SL.earth","SL.gam",
                                                     "SL.glm","SL.glmnet",
                                                     "SL.glm.interaction",
                                                     "SL.mean")){
  # following Edward's implementation of semi-parametric density estimation.
  pimod <- SuperLearner(Y = A, X = data.frame(W), SL.library = SL.library)
  pimod.vals <- pimod$SL.predict
  pi2mod <- SuperLearner(Y = (A - pimod.vals)^2, X = data.frame(W),
                         SL.library = SL.library)
  pi2mod.vals <- pi2mod$SL.predict
  a.std <- (A - pimod.vals)/sqrt(pi2mod.vals)
  get.pi <- function(a0, w0){
    pred.mean <- predict(pimod, data.frame(w0))$pred
    pred.var <- predict(pi2mod, data.frame(w0))$pred
    a.std.pred <- (a0 - pred.mean)/sqrt(pred.var)
    dnorm(a.std.pred)
  }
  a.std <- (A - mean(A))/sd(A)
  density.est <- density(a.std)
  marginal.density <- approxfun(density.est$x, density.est$y)
  get.marginal <- function(a){
    a.std.pred <- (a - mean(A))/sd(A)
    dnorm(a.std.pred)
  }
  function(a0, w0) get.pi(a0, w0)/get.marginal(a0)
}

#' Title: .fit.density
#'
#' Estimate standardized propensity score using SuperLearner-based model.
#' The implementation is by Ted and the original method was proposed in the
#' paper by Diaz and van der Laan
#'
#' @param A \code{n x 1} numeric vector of observed exposure values.
#' @param W \code{n x p} matrix of potential confounders.
#' @param SL.library list of base learners for SuperLearner package
#' @param cdf TRUE if it transforms the exposure A to \code{[0,1]} by empirical CDF
#'  the default is TRUE.
#' @export
#' @return
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
                       n.bins=2:floor(length(U)/50))
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
