# Per-subject prediction and Jacobian for mixed-effects engines

Unlike the population gradient solver, which applies a single `theta` to
every subject, this takes an `nsub x ntheta` matrix whose row `id` holds
that subject's parameter vector (`phi = beta + b`, as supplied by
`lme4::nlmer`). Each subject is solved reusing the loaded `thetaGrad`
model, sticky-tolerance recalculation, event finite differences, and
jump sensitivities. The nonlinear problem must already be loaded with
[`.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.md).

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
