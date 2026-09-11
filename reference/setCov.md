# Set the covariance type based on prior calculated covariances

Switches a completed fit's covariance to `method`. A previously computed
covariance is re-installed from the cache; otherwise it is recomputed at
the converged estimates: `"r,s"`/`"r"`/`"s"` and `"analytic"` on a
zero-iteration FOCEI model, and `"sa"` (SAEM Louis FIM) / `"imp"`
(importance-sampling Monte-Carlo) via the decoupled recompute engine
(the latter two require a mixed-effects fit). When
`"sa"`/`"imp"`/`"analytic"` cannot be computed the covariance is left
unchanged (it is never silently downgraded to `"r,s"`).

## Usage

``` r
setCov(fit, method)
```

## Arguments

- fit:

  nlmixr2 fit

- method:

  covariance method (see the \`covMethod\` argument for the control
  options for the choices)

## Value

Fit object with covariance updated

## Details

Every focei covariance comes in two shapes (see `covFull` in
[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md)),
and both are named: `"r,s"`, `"r"`, `"s"` and `"analytic"` are the
structural-theta block, while `"r,s (full)"`, `"r (full)"`, `"s (full)"`
and `"analytic (full)"` are the full theta + residual sigma + Omega
matrix. A focei fit computes both and caches the one it does not
install, so swapping between them costs nothing. The shapes are not
submatrices of one another on the finite-difference path – `"s"` inverts
the theta block of the cross-product while `"s (full)"` takes the theta
block of the full inverse, which also carries the Omega estimation
uncertainty – so the standard errors differ. On the analytic path the
assembly is always full and `"analytic"` is a submatrix of
`"analytic (full)"`, so the theta standard errors agree.

`fit$covMethod` names the installed covariance and `names(fit$covList)`
the cached alternatives (the fit print shows both).

## See also

[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Author

Matt Fidler
