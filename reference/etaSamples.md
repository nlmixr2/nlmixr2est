# Draw each subject's etas from their conditional distribution

Draws \`n\` samples per subject from the conditional (posterior)
distribution of the random effects, \`p(eta \| y, theta, Omega,
sigma)\`, at the fit's final population estimates. Unlike the EBEs, the
draws carry each subject's uncertainty, so they can be used to simulate
that subject rather than a new one.

## Usage

``` r
etaSamples(fit, n = 500L, method = c("auto", "laplace", "imp"), seed = 42L)
```

## Arguments

- fit:

  nlmixr2 fit with random effects

- n:

  number of draws per subject

- method:

  sampler: \`"auto"\`, \`"laplace"\` or \`"imp"\`

- seed:

  random seed; the same seed gives the same draws whatever the number of
  cores

## Value

a data frame with columns \`ID\`, \`draw\` and one column per eta,
ordered by \`draw\` then \`ID\`. Passed as \`params\` to
\`rxode2::rxSolve(fit, params = draws, nStud = n, simVariability =
FALSE)\` it simulates each subject once per draw (\`sim.id\` is the
draw). Its \`"method"\` attribute names the sampler.

## Details

\`method = "laplace"\` draws from a normal centered at each subject's
conditional mode with the inverse of the inner Hessian (\`fit\$phiC\`)
as its covariance. A fit whose etas are not conditional modes (for
example SAEM or emvi) gets both from a posthoc FOCEi step at its
estimates. \`method = "imp"\` resamples the importance samples of
\[addCondDist()\]. \`method = "auto"\` is \`"laplace"\`.

Mixture, IOV and non-parametric fits are not supported, nor is an Omega
that is not positive definite.

## See also

\[addCondDist()\]

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
#> → Calculating residuals/tables
#> ✔ done
smp <- etaSamples(fit, n = 100)
# each subject's conditional SD of eta.cl
aggregate(eta.cl ~ ID, smp, sd)
#>    ID     eta.cl
#> 1   1 0.10299118
#> 2   2 0.10098505
#> 3   3 0.10004002
#> 4   4 0.10274447
#> 5   5 0.09410539
#> 6   6 0.11230115
#> 7   7 0.10968900
#> 8   8 0.08755851
#> 9   9 0.09963759
#> 10 10 0.11082492
#> 11 11 0.10010616
#> 12 12 0.09257853
# simulate each subject with its own uncertainty
s <- rxode2::rxSolve(fit, params = smp, nStud = 100, simVariability = FALSE)
#>  
#>  
#> ℹ using original fit data for simulation
#> Warning: multi-subject simulation without 'omega'
# }
```
