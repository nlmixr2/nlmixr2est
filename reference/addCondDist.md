# Add each subject's conditional distribution of the etas to a fit

Samples each subject's conditional (posterior) distribution of the
random effects, \`p(eta \| y, theta, Omega, sigma)\`, at the fit's final
population estimates, and stores its mean and SD as \`fit\$condDist\`.
It applies to a fit from any mixed-effects estimation method.

## Usage

``` r
addCondDist(
  fit,
  method = c("auto", "imp", "saem"),
  nsamp = 300L,
  keepSamples = FALSE,
  burn = 50L,
  seed = 42L
)
```

## Arguments

- fit:

  nlmixr2 fit with random effects

- method:

  the sampler: \`"auto"\` (default) is \`"saem"\` for a SAEM fit and
  \`"imp"\` otherwise

- nsamp:

  number of samples per subject (default 300)

- keepSamples:

  Boolean (default \`FALSE\`); keep \`nsamp\` equally weighted draws per
  subject as the \`"samples"\` attribute of \`fit\$condDist\`, an \`ID x
  eta x draw\` array

- burn:

  MCMC iterations discarded by a held \`"saem"\` run (default 50)

- seed:

  random seed of the sampler

## Value

the fit, with \`fit\$condDistShrinkage\`, the eta shrinkage ( the
conditional variances, \`100 \* mean(SD^2) / omega\` (an unsampled
subject's \`SD^2\` taken as \`omega\`), and \`fit\$condDist\`: a data
frame \`ID\`, \`ETA\`, \`MEAN\`, \`SD\`. Its \`"method"\` attribute
names the sampler (\`"imp"\`, \`"saem"\` for a SAEM fit's own chains or
\`"mcmc"\` for a held SAEM run) and, for \`"imp"\`, its \`"ess"\`
attribute holds each subject's effective sample size.

## Details

\`method = "imp"\` draws importance samples from a proposal centered at
each subject's mode, holding every population parameter at the fit's
estimate (the frozen \`est = "imp"\` E-step), and weights them by the
posterior. Mixture, IOV, non-parametric and repeated (\`same()\`) Omega
block fits are not supported.

\`method = "saem"\` uses the SAEM MCMC kernels. For a SAEM fit, the
chains of its own last iterations are the draws when they are enough
(\`nsamp\` draws from \`nmc\` chains, within the last half of the
second-phase iterations, where the parameters have nearly settled), so
nothing is simulated; otherwise the kernels run with every population
parameter held at the fit's estimates, continuing a SAEM fit's chains
or, for any other fit, starting cold, and \`burn\` iterations are
discarded; it needs a model SAEM can fit (every eta on a population
parameter of its own). An eta SAEM does not sample for a subject
(uninformative) has \`SD\` \`NA\` and \`MEAN\` 0, its prior, as in
\`fit\$eta\`; that is exact when the eta is uncorrelated with the
others. A model SAEM fits with temporary etas (a likelihood parameter
with no eta) is sampled with \`"imp"\`, since those chains also sampled
parameters the reported model holds fixed.

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
