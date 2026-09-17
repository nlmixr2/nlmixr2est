# Set the covariance type based on prior calculated covariances

Switches a completed fit's covariance to `method`. A previously computed
covariance is re-installed from the cache; otherwise it is recomputed at
the converged estimates: `"r,s"`/`"r"`/`"s"` and `"analytic"` on a
zero-iteration FOCEI model, and `"sa"` (SAEM Louis FIM) / `"imp"`
(importance-sampling Monte-Carlo) via the decoupled recompute engine
(the latter two require a mixed-effects fit). When a covariance cannot
be computed it is left unchanged (it is never silently downgraded to
`"r,s"`).

## Usage

``` r
setCov(fit, method, ...)

# Default S3 method
setCov(fit, method, ...)

# S3 method for class 'analytic'
setCov(fit, method, ...)

# S3 method for class '`r,s`'
setCov(fit, method, control = rsControl(), ...)

# S3 method for class 'r'
setCov(fit, method, control = rsControl(), ...)

# S3 method for class 's'
setCov(fit, method, control = rsControl(), ...)

# S3 method for class 'sa'
setCov(fit, method, control = saControl(), ...)

# S3 method for class 'imp'
setCov(fit, method, control = impCovControl(), ...)
```

## Arguments

- fit:

  nlmixr2 fit

- method:

  covariance method (see the `covMethod` argument of the control options
  for the choices)

- ...:

  arguments passed to the covariance method

- control:

  options for the covariance method itself, only needed to change its
  defaults:
  [`rsControl()`](https://nlmixr2.github.io/nlmixr2est/reference/rsControl.md)
  for `"r,s"`, `"r"` and `"s"`,
  [`saControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saControl.md)
  for `"sa"` and
  [`impCovControl()`](https://nlmixr2.github.io/nlmixr2est/reference/impCovControl.md)
  for `"imp"`

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

## Adding a covariance method

`setCov()` is an S3 generic dispatched on `method` (without any
`" (full)"` suffix), so another package (for example SIR or a bootstrap)
can add a covariance by registering a method:


    # registered with S3method(nlmixr2est::setCov, sir) in NAMESPACE
    setCov.sir <- function(fit, method, control = sirControl(), ...) {
      # compute on the estimation scale, named like fit$cov
      mySirCovariance(fit, control)
    }

A method that has options declares its own `control` holding only those
options, so `setCov(fit, "sir")` uses the defaults and
`setCov(fit, "sir", control = sirControl(...))` changes them. The method
receives the fit and `method` (carrying dispatch classes; use
`unclass(method)` for the plain name) and either returns a named
covariance matrix, which `setCov()` checks for positive definiteness and
installs as `method` (updating the standard errors and keeping the prior
covariance in `fit$covList`), or installs the covariance itself and
returns `NULL`. A mixture fit's matrix is rotated onto the probability
scale unless it carries `attr(, "mixRotated")` set to `TRUE`. A method
that cannot compute the covariance should
[`stop()`](https://rdrr.io/r/base/stop.html).
[`setCovAllMethods()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovAllMethods.md)
lists the available methods.

Each covariance remembers the options it was computed with (in
`fit$env$covOptions`). A covariance already on the fit is reinstalled
from `fit$covList` only when the requested options – the supplied
`control`, or the method's default one – are the same; otherwise it is
recomputed. A covariance computed during estimation used the fit's own
settings.

## See also

[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Author

Matt Fidler
