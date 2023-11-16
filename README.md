# DebiasedDoseResponse

Code repository for the paper "Debiased inference for a covariate-adjusted regression function". See <https://arxiv.org/abs/2210.06448> for details.

The numerical studies in the original paper can be found [here](https://github.com/Kenta426/sim-debiased-inference).

To install this `R` package, first install the `devtools` package. Then type:

    library(devtools)
    devtools::install_github("Kenta426/DebiasedDoseResponse")
    library(DebiasedDoseResponse)

## Usage

```r
set.seed(10000)
library(DebiasedDoseResponse)
library(ggplot2)

n <- 500; cols <- 3
W <- matrix(runif(n*cols), ncol = cols)   # a 500 * 3 matrix of covariates
A <- rnorm(n, mean=W%*%rnorm(cols))       # a 500 * 1 vector of treatment variable
Y <- rnorm(n, mean = sin(2*A)+W[,1])      # a 500 * 1 vector of response variable
est.res <- debiased_inference(Y, A, W)    # compute debiased local linear 
p <- plot_debiased_curve(est.res)         # plot debiased local linear 

# compute true covariate-adjusted regression by integrating W
g <- function(a){
  points <- function(a0){
    integrate(function(t){sin(2*a0)+t}, lower=0, upper=1)$value}
  sapply(a, points)
}
# plot true covariate-adjusted regression in coral
p + geom_line(aes(x=eval.pts, y=g(eval.pts)), color="coral")
```

<p align="center">
  <img src="https://github.com/Kenta426/DebiasedDoseResponse/blob/main/figs/demo1.png" />
</p>


```r
Y <- rnorm(n, mean = sin(2*A*W[,1])) 

est.res <- debiased_inference(Y, A, W) 
p <- plot_debiased_curve(est.res)

# compute true covariate-adjusted regression by integrating W
g <- function(a){
  points <- function(a0){
    integrate(function(t){sin(2*a0*t)}, lower=0, upper=1)$value}
  sapply(a, points)
}
p + geom_line(aes(x=eval.pts, y=g(eval.pts)), color="coral")
```
<p align="center">
  <img src="https://github.com/Kenta426/DebiasedDoseResponse/blob/main/figs/demo2.png" />
</p>

## Recommendation
For large data, it is recommended to compute outcome regression and standardized
density separately. The example is provided below. This avoids computing nuisance 
estimators every time.
```r
set.seed(10000)
library(DebiasedDoseResponse)
library(ggplot2)

n <- 500; cols <- 3
W <- matrix(runif(n*cols), ncol = cols)   # a 500 * 3 matrix of covariates
A <- rnorm(n, mean=W%*%rnorm(cols))       # a 500 * 1 vector of treatment variable
Y <- rnorm(n, mean = sin(2*A)+W[,1])      # a 500 * 1 vector of response variable

mu.hat <- .fit.regression(Y, A, W)         # fit outcome regression 
g.hat <- .fit.semiparametric.density(A, W) # fit standardized propensity  

est.res <- debiased_inference(Y, A, W, mu=mu.hat, g=g.hat)    # compute debiased local linear 
```

