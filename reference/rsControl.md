# Options for the finite-difference covariance in setCov()

Used by `setCov(fit, "r,s")`, `"r"` and `"s"`. Each option left `NULL`
keeps the value the fit was estimated with.

## Usage

``` r
rsControl(
  hessEps = NULL,
  gillKcov = NULL,
  gillStepCov = NULL,
  gillFtolCov = NULL,
  covGillF = NULL,
  covSmall = NULL,
  rmatNorm = NULL,
  smatNorm = NULL
)
```

## Arguments

- hessEps:

  is a double value representing the epsilon for the Hessian
  calculation. This is used for the R matrix calculation.

- gillKcov:

  Max steps to determine the optimal forward/central difference step
  size per parameter (Gill 1983) during the covariance step. \`0\` = no
  optimal step size determined.

- gillStepCov:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration during the covariance step is equal to the new step size =
  (prior step size)\*gillStepCov

- gillFtolCov:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates during the
  covariance step.

- covGillF:

  Use the Gill calculated optimal Forward difference step size for the
  instead of the central difference step size during the central
  difference gradient calculation.

- covSmall:

  Small number used to compare covariance estimates (sandwich vs R/S
  matrix) before rejecting one as too small to be the final covariance
  estimate.

- rmatNorm:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the R matrix

- smatNorm:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the S matrix

## Value

`rsControl` object

## See also

[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md)

## Author

Matt Fidler

## Examples

``` r
rsControl(hessEps = 1e-4)
#> $hessEps
#> [1] 1e-04
#> 
#> attr(,"class")
#> [1] "rsControl"
```
