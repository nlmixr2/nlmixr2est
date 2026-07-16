# Recover and finalize the resident nlm parameter history

Returns the parameter history accumulated in the resident nlm scaling
struct (one row per iteration type per recorded evaluation) as a data
frame, and stops further recording/printing (`save` and `every` are
reset to 0). Must be called while
[`.nlmSetupEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmSetupEnv.md)
is still loaded – i.e. before
[`.nlmFreeEnv()`](https://nlmixr2.github.io/nlmixr2est/reference/dot-nlmFreeEnv.md).
Used by `.nlmFinalizeList` for the standard nlm-family estimators and
directly by externally-optimized engines such as `babelmixr2`'s nlmer.

## Usage

``` r
nlmGetParHist(p = TRUE)
```

## Arguments

- p:

  When `TRUE` (default) also print the final iteration line.

## Value

A data frame of the recorded parameter history.

## Details

This is an internal function and should not be called directly.

## Author

Matthew L. Fidler
