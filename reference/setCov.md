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

## See also

[`foceiControl()`](https://nlmixr2.github.io/nlmixr2est/reference/foceiControl.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Author

Matt Fidler
