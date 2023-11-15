###############################################################################
# Super-learning a conditional density
#
# The functions in this file implement super-learning of a conditional density
# based on a user-specified library of estimators. The methods are based in part on:
#
# D?az Munoz, Ivan, and Mark J. van der Laan. "Super learner based conditional density estimation with application to marginal structural models." The international journal of biostatistics 7.1 (2011): 1-20
#
# Author: tedwestling
###############################################################################

########################################################
# Function: CD.SuperLearner
#
# This function performs Super-Learning estimation
# of a conditional density p(x | w), where x is a continuous
# scalar response and w is a possibly multivariate set of
# covariates.
#
# Arguments:
# X - a numeric vector of observations of the variable
#     whose conditional density is to be estimated
# W - a data frame of numeric conditioning covariates
# SL.library - the library of candiate algorithms to use.
#     This should have the same format as required by the
#     SuperLearner package, and any of the algorithms in that
#     package can be specified here as well. This library will
#     be used when estimating the conditional mean and variance
#     functions and when estimating the bin probabilities (see
#     details).
# n.folds - the number of folds for the cross-validation. Defaults
#     to 10
# n.bins - an integer vector of the number of bins to use for
#     the discretized procedures.
# verbose - an indicator of whether the function should output
#     progress to the command line
#
# Value:

# Details: This function uses a meta-super learner to combine two
#     types of conditional density estimation. The first method
#     is a "location-scale" method also used in the csteff function
#     in the npcausal package. First the conditional mean of X given
#     W is estimated using SuperLearner, then the SL fitted conditional
#     mean is used to estimate the conditional variance of X given W using
#     a second SuperLearner -- i.e. the conditinoal mean of
#     (X - fitted conditional mean)^2. Finally, a univariate kernel smoother
#     is applied to the standardized values
#     (X - fitted conditional mean) / sqrt(fitted conditional variance).
#     This method works well when the conditional density is approximately
#     p(x | w) = f((x - mu(x | w)) / sigma(x | w)) for some univariate f.
#     The second method is more nonparametric method based on binning. It is
#     an adaptation of the algorithm described in:
#     D?az Munoz, Ivan, and Mark J. van der Laan. "Super learner based conditional density estimation with application to marginal structural models.
#     For a given number of bins k, x is discretized in to k equal-weight bins.
#     The conditional probability of x being in each bin given w is estimated
#     using the supplied super learning library, and these probabilities are then
#     normalized to sum to one. The conditional density p(x | w) is then approximated by
#     the conditional probability of the bin that x falls in divided by the width of
#     the bin.
#     A "meta" super learner is used to combine the location-scale
#     method and all of the bin methods to a single estimate.

CD.SuperLearner <- function(X, W, SL.library, n.folds=10,
                            n.bins=2:floor(length(X)/50),
                            verbose=FALSE, save.threshold=.001) {
  n <- length(X)
  sorted  <- rep(1:n.folds, length.out = n)
  folds <- sample(sorted, n, replace=FALSE)
  valid.rows <- lapply(1:n.folds, function(v) which(folds == v))

  # Construct prediction matrix
  # preds will contain the test-data predictions from each model on each bin type.
  if(class(SL.library) == "list") {
    cand.algs <- unlist(lapply(SL.library, function(vec) paste(vec[1], vec[-1], sep='_')))
  } else {
    cand.algs <- paste(SL.library, "All", sep="_")
  }
  n.sl.models <- length(cand.algs)
  n.models <- 1 + length(n.bins) * n.sl.models
  library.densities <- cv.library.densities <- matrix(NA, nrow=n, ncol=n.models)
  library.names <- c("SL.loc.scale", c(sapply(n.bins, function(k) paste0(cand.algs, "_", k, "bins"))))

  # model_fits will contain each of the models fitted to the entire dataset
  model.fits <- vector(mode='list', 1+length(n.bins))

  # First fit the nonparametric "location-scale" smoothing-based model
  if(verbose) cat("Fitting location-scale... ")
  for(v in 1:n.folds) {
    if(verbose) cat("fold", v, "... ")
    capture.output(mean.model <- try(SuperLearner(Y=X[folds != v], X=W[folds != v,],
                                                  SL.library = SL.library,
                                                  family='gaussian',
                                                  method="method.NNLS"), silent=TRUE))
    if(class(mean.model) != "try-error") {
      mean.preds <- predict(mean.model, newdata=W[folds == v,], onlySL = TRUE)$pred
      capture.output(var.model <- try(SuperLearner(Y=(X[folds != v] - mean.model$SL.predict)^2,
                                                   X=W[folds != v,],
                                                   SL.library = SL.library,
                                                   family='gaussian',
                                                   method="method.NNLS"), silent=TRUE))
      if(class(var.model) != "try-error") {
        var.preds <- predict(var.model, newdata=W[folds == v,], onlySL = TRUE)$pred
        var.preds[var.preds <= 0] <- min(var.preds[var.preds > 0])
        var.hats <- c(var.model$SL.predict )
        var.hats[var.hats <= 0] <- min(var.hats[var.hats > 0])
        if(any(is.na((X[folds != v] - mean.model$SL.predict)/sqrt(var.hats)))) {
          stop("missing values error")
        }
        dens <- density((X[folds != v] - mean.model$SL.predict)/sqrt(var.hats))
        cv.library.densities[folds == v, 1] <- c(approx(dens$x, dens$y, xout=(X[folds == v] - mean.preds)/sqrt(var.preds), rule=2)$y / sqrt(var.preds))
      }

    }
  }

  if(verbose) cat("full data... ")
  mean.model <- SuperLearner(Y=X, X=W, SL.library = SL.library,
                             family='gaussian', method="method.NNLS")
  mean.preds <- mean.model$SL.predict
  var.model <- SuperLearner(Y=(X - mean.preds)^2, X=W, SL.library = SL.library,
                            family='gaussian', method="method.NNLS")
  var.preds <- var.model$SL.predict
  var.preds <- pmax(var.preds, 1e-5)
  dens <- density((X - mean.preds)/sqrt(var.preds))
  library.densities[, 1] <- c(approx(dens$x, dens$y, xout=(X - mean.preds)/sqrt(var.preds))$y / sqrt(var.preds))
  model.fits[[1]] <- list(mean.model=mean.model, var.model=var.model, dens=dens)

  # column of prediciton matrix counter
  column.start <- 2
  column.end <- 1 + n.sl.models
  bin.breaks <- NULL
  # Each k defines a number of equal-area bins for the histogram
  for(j in 1:length(n.bins)) {
    k <- n.bins[j]
    if(verbose) cat("\nEstimating models with", k, "bins... ")
    ycuts <- seq(0,1,by=1/k)
    bins <- quantile(X, ycuts)

    disc.X <- cut(X, breaks=bins, include.lowest=TRUE, ordered_result=TRUE)
    disc.X.num <- as.numeric(disc.X)
    labs <- levels(disc.X)
    cuts <- cbind(bins[-length(bins)], bins[-1])
    bin.lengths <- as.numeric(cuts[,2] - cuts[,1])
    bin.breaks[[paste0("n.bins.", k)]] <- list(breaks=bins, bin.lengths=bin.lengths)
    cv.bin.probs <- bin.probs <- array(NA, dim=c(n, n.sl.models, k))
    for(bin in 1:k) {
      if(verbose) cat("bin", bin, "... ")
      capture.output(bin.fit <- try(SuperLearner(Y=as.numeric(disc.X.num==bin), X=W, family='binomial', SL.library = SL.library, method='method.NNloglik', cvControl = list(V=n.folds, validRows=valid.rows)), silent=TRUE))
      if(class(bin.fit) != "try-error") {
        model.fits[[1+j]][[paste0("bin", bin, ".model")]] <- bin.fit
        cv.bin.probs[,,bin] <- bin.fit$Z
        bin.probs[,,bin] <- bin.fit$library.predict
      }
    }
    for(i in 1:n) {
      for(m in 1:n.sl.models) {
        cv.bin.probs[i,m,] <- cv.bin.probs[i,m,] / sum(cv.bin.probs[i,m,] )
        bin.probs[i,m,] <- bin.probs[i,m,] / sum(bin.probs[i,m,] )
      }
      cv.library.densities[i,column.start:column.end] <- cv.bin.probs[i,,disc.X.num[i]] / bin.lengths[disc.X.num[i]]
      library.densities[i,column.start:column.end] <- bin.probs[i,,disc.X.num[i]] / bin.lengths[disc.X.num[i]]
    }

    column.start <- column.end + 1
    column.end <- column.end + n.sl.models
  }

  if(verbose) cat("\nOptimizing model weights...\n")

  # Remove algs with errors in cv predictions
  errors.in.library <- apply(cv.library.densities, 2, function(col) any(is.na(col)))
  if(any(errors.in.library)) warning(paste0("Errors in the following candidate algorithms: ", library.names[which(errors.in.library)]))
  n.include <- sum(!errors.in.library)

  # Do SL log-likelihood optimization
  cv_risk <- function(beta) -sum(log(cv.library.densities[,!errors.in.library] %*% beta))
  capture.output(solnp_solution <- solnp(rep(1/n.include, n.include), cv_risk, eqfun=sum, eqB=1, ineqfun=function(beta) beta, ineqLB=rep(0,n.include), ineqUB=rep(1, n.include)))
  coef <- rep(0, n.models)
  coef[!errors.in.library] <- solnp_solution$pars
  SL.density <- c(library.densities[,!errors.in.library] %*% solnp_solution$pars)

  # Only return algorithms with larger than save.threshold coefficient
  save <- which(coef > save.threshold)
  save.model.fits <- vector(mode='list', n.models)
  model.bins <- c(NA, rep(n.bins, each=n.sl.models))
  model.algs <- c(NA, rep(cand.algs, length(n.bins)))
  for(j in save) {
    if(j == 1) save.model.fits[[1]] <- model.fits[[1]]
    else {
      this.bins <- model.bins[j]
      this.alg <- model.algs[j]
      for(bin in 1:this.bins) {
        save.model.fits[[j]][[paste0("bin", bin, ".model")]] <- model.fits[[1 + which(n.bins == this.bins)]][[paste0("bin", bin, ".model")]]$fitLibrary[[this.alg]]
      }
    }
  }

  return(list(SL.library=SL.library, library.names=library.names, model.bins=model.bins, model.algs=model.algs, SL.density=SL.density, library.densities=library.densities, coef=coef, cv.library.densities=cv.library.densities, model.fits=save.model.fits, n.bins=n.bins, bin.breaks=bin.breaks, n.folds=n.folds, folds=folds, errors.in.library=errors.in.library))
}

predict.CD.SuperLearner <- function(fit, newdata, threshold=.001, verbose=FALSE) {
  model.nums <- which(fit$coef > threshold & !unlist(lapply(fit$model.fits, is.null)))
  new.coef <- fit$coef[model.nums]
  new.coef <- new.coef / sum(new.coef)
  model.names <- fit$library.names[model.nums]
  bins <- fit$model.bins[model.nums]
  bins <- unique(bins[!is.na(bins)])

  newW <- newdata[,names(newdata) != "X"]
  library.densities <- matrix(NA, nrow=nrow(newdata), ncol=length(fit$library.names))

  for(k in model.nums) {
    if(k == 1) {
      mean.preds <- c(predict(fit$model.fits[[1]]$mean.model, newdata=newW, onlySL=TRUE)$pred)
      var.preds <- c(predict(fit$model.fits[[1]]$var.model, newdata=newW, onlySL=TRUE)$pred)
      library.densities[,1] <- c(approx(fit$model.fits[[1]]$dens$x, fit$model.fits[[1]]$dens$y, xout=(newdata$X - mean.preds)/sqrt(var.preds))$y / sqrt(var.preds))
    } else {
      bins <- fit$model.bins[k]
      alg <- fit$model.algs[k]
      disc.X <- findInterval(newdata$X, fit$bin.breaks[[paste0("n.bins.", bins)]]$breaks, all.inside = TRUE)
      bin.probs <- matrix(NA, nrow=nrow(newdata), ncol=bins)
      for(bin in 1:bins) {
        bin.probs.bin <- try(predict(fit$model.fits[[k]][[paste0("bin", bin, ".model")]], newdata=newW, family='binomial'), silent=TRUE)
        if(class(bin.probs.bin) != "try-error") {
          bin.probs[,bin] <- bin.probs.bin
        } else {
          bin.probs[,bin] <- NA
        }
      }
      for(i in 1:nrow(newdata)) {
        bin.probs[i,] <- bin.probs[i,] / sum(bin.probs[i,] )
        library.densities[i,k] <- bin.probs[i,disc.X[i]] / fit$bin.breaks[[paste0("n.bins.", bins)]]$bin.lengths[disc.X[i]]
      }
    }
  }
  col.has.na <- apply(library.densities[,model.nums, drop=FALSE], 2, function(col) any(is.na(col)))
  if(sum(!col.has.na) == 0) stop("All algorithms had errors in prediction. Please try a different library.")
  new.coef <- new.coef[!col.has.na] / sum(new.coef[!col.has.na])
  list(library.densities=library.densities, SL.density=c(library.densities[,model.nums[!col.has.na], drop=FALSE] %*% new.coef))
}
