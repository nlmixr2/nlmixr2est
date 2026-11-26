# Add each subject's conditional distribution of the etas to a fit

Samples each subject's conditional (posterior) distribution of the
random effects, \`p(eta \| y, theta, Omega, sigma)\`, at the fit's final
population estimates, and stores its mean and SD as \`fit\$condDist\`.
It applies to a fit from any mixed-effects estimation method.

## Usage

``` r
addCondDist(
  fit,
  method = c("auto", "imp"),
  nsamp = 300L,
  keepSamples = FALSE,
  seed = 42L
)
```

## Arguments

- fit:

  nlmixr2 fit with random effects

- method:

  \`"auto"\` (default) or \`"imp"\`: the sampler

- nsamp:

  number of samples per subject (default 300)

- keepSamples:

  Boolean (default \`FALSE\`); keep \`nsamp\` equally weighted draws per
  subject as the \`"samples"\` attribute of \`fit\$condDist\`, an \`ID x
  eta x draw\` array

- seed:

  random seed of the sampler

## Value

the fit, with \`fit\$condDist\`: a data frame \`ID\`, \`ETA\`, \`MEAN\`,
\`SD\`. Its \`"method"\` attribute names the sampler and its \`"ess"\`
attribute holds each subject's effective sample size.

## Details

\`method = "imp"\` draws importance samples from a proposal centered at
each subject's mode, holding every population parameter at the fit's
estimate (the frozen \`est = "imp"\` E-step), and weights them by the
posterior. Mixture, IOV, non-parametric and repeated (\`same()\`) Omega
block fits are not supported.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{
one.cmt <- function() {
  ini({
    tka <- 0.45
    tcl <- 1
    tv <- 3.45
    eta.ka ~ 0.6
    eta.cl ~ 0.3
    eta.v ~ 0.1
    add.sd <- 0.7
  })
  model({
    ka <- exp(tka + eta.ka)
    cl <- exp(tcl + eta.cl)
    v <- exp(tv + eta.v)
    linCmt() ~ add(add.sd)
  })
}
fit <- nlmixr2(one.cmt, nlmixr2data::theo_sd, est = "focei", control = list(print = 0))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → calculate ∂(f)/∂(η)
#> → calculate ∂(R²)/∂(η)
#> → finding duplicate expressions in inner model...
#> → finding duplicate expressions in EBE model...
#> → compiling inner model...
#>  
#>  
#> ✔ done
#> → finding duplicate expressions in FD model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → compiling events FD model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
fit <- addCondDist(fit)
head(fit$condDist)
#>   ID    ETA         MEAN         SD
#> 1  1 eta.ka  0.086724343 0.13290765
#> 2  1 eta.cl -0.485157088 0.09737933
#> 3  1  eta.v -0.089112392 0.04531539
#> 4  2 eta.ka  0.198230689 0.15866379
#> 5  2 eta.cl  0.144418922 0.09494652
#> 6  2  eta.v  0.005335317 0.06765941
# }
```
