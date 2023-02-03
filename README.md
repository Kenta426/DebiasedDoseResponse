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
n <- 500; cols <- 3
W <- matrix(runif(n*cols), ncol = cols) # a 200 * 3 matrix of covariates
A <- rnorm(n, mean=W%*%rnorm(cols)) # a 200 * 1 vector of treatment variable
Y <- rnorm(n, mean = sin(2*A)+W[,1]) # a 200 * 1 vector of response variable
est.res <- debiased_inference(Y, A, W)
p <- plot_debiased_curve(est.res)
p + geom_line(aes(x=eval.pts, y=sin(eval.pts*2)+1/2), color="coral")
```

<p align="center">
  <img src="https://github.com/Kenta426/DebiasedDoseResponse/blob/main/figs/demo1.png" />
</p>
