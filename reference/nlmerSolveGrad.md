# Per-subject prediction and Jacobian for mixed-effects engines

Like the population gradient solver but takes a per-subject
`nsub x ntheta` parameter matrix (`phi = beta + b`, as supplied by
`lme4::nlmer`) instead of one shared `theta`. Requires
[`.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.md)
to already be loaded.

## Usage

``` r
nlmerSolveGrad(thetaMat)
```

## Arguments

- thetaMat:

  A `nsub x ntheta` matrix of per-subject parameter values. Row `id` is
  solved against subject `id` (in the loaded `etTrans` order).

## Value

A `nobsTot x (ntheta+1)` matrix in the loaded (`etTrans`) observation
order: column 1 is the prediction (`rx_pred_`) and columns 2..(ntheta+1)
are `d(pred)/d(THETA[i])`.

## Details

This is an internal function and should not be called directly.

## Author

Matthew L. Fidler
