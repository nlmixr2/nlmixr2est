# Unscale an nlm-family parameter vector back to the natural scale

Converts a parameter vector from the estimation (scaled) scale used by
the currently loaded nlm-family problem back to the natural scale, using
the scaling that `nlmSetup()` installed. Exported so that external
engines driving the nlm-family objective (e.g. `babelmixr2`'s FME-based
methods) do not have to reach into the namespace for it (#940).

## Usage

``` r
nlmUnscalePar(p)
```

## Arguments

- p:

  Numeric parameter vector on the estimation scale; its length must
  match the loaded problem's number of parameters.

## Value

A numeric vector of the same length (names preserved) on the natural
scale.

## Details

An nlm-family problem must be loaded (via the internal
[`.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.md)/`nlmSetup()`
path) when this is called; the scaling is part of that problem's state.

## Author

Matthew L. Fidler
