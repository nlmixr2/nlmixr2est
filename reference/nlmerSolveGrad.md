# Per-subject prediction and Jacobian for mixed-effects engines

Like the population gradient solver but takes a per-subject
`nsub x ntheta` parameter matrix (`phi = beta + b`, as supplied by
`lme4::nlmer`) instead of one shared `theta`. Requires
[`.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.md)
to already be loaded.

## Usage

``` r
nlmerSolveGrad(thetaMat, record = FALSE)
```

## Arguments

- thetaMat:

  A `nsub x ntheta` matrix of per-subject parameter values. Row `id` is
  solved against subject `id` (in the loaded `etTrans` order).

- record:

  When `TRUE`, record this evaluation's population parameter estimate –
  the per-subject mean of `thetaMat`'s columns (`phi = beta + b`
  averaged over subjects, which equals the fixed effect exactly for
  parameters without a random effect) – into the resident nlm parameter
  history via the shared scale machinery. This is how an external
  optimizer (e.g. `lme4::nlmer`) populates the iteration print and the
  history recovered by
  [`nlmGetParHist()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmGetParHist.md).
  No objective value is recorded (the scale's `showOfv` is expected to
  be 0 for these engines). Defaults to `FALSE`.

## Value

A `nobsTot x (ntheta+1)` matrix in the loaded (`etTrans`) observation
order: column 1 is the prediction (`rx_pred_`) and columns 2..(ntheta+1)
are `d(pred)/d(THETA[i])`.

## Details

This is an internal function and should not be called directly.

## Author

Matthew L. Fidler
