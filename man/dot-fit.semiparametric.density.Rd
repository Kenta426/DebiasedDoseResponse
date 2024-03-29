% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_nuisances.R
\name{.fit.semiparametric.density}
\alias{.fit.semiparametric.density}
\title{Title: .fit.semiparametric.density
Estimate standardized propensity score using semiparametric model.
The implementation is based on npcausal package by Edward Kennedy.}
\usage{
.fit.semiparametric.density(
  A,
  W,
  SL.library = c("SL.earth", "SL.gam", "SL.glm", "SL.glmnet", "SL.glm.interaction",
    "SL.mean")
)
}
\arguments{
\item{A}{\code{n x 1} numeric vector of observed exposure values.}

\item{W}{\code{n x p} matrix of potential confounders.}

\item{SL.library}{list of base learners for SuperLearner package}
}
\value{

}
\description{
Title: .fit.semiparametric.density
Estimate standardized propensity score using semiparametric model.
The implementation is based on npcausal package by Edward Kennedy.
}
\examples{
# Sample data
n <- 200; cols <- 3
W <- matrix(runif(n*cols), ncol=cols) # 200 * 3 matrix of covariates
A <- rnorm(n, mean=W\%*\%rnorm(cols))   # 200 * 1 vector of treatment variable
Y <- rnorm(n, mean=sin(A)+W\%*\%rnorm(cols))  # 200 * 1 vector of response
mu.hat <- .fit.semiparametric.density(A, W)
mu.hat(A, W) # estimated values
}
