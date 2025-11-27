# This is an internal function for modifying the UI to apply mu2 referencing

mu2 referencing is algebraic mu-referencing by converting to the
transformation to a single value in the original dataset, and moving
that around

## Usage

``` r
.uiApplyMu2(env)
```

## Arguments

- env:

  Environment needed for nlmixr2 fits

## Value

Either the original model() block (if changed) or NULL if not changed

## Author

Matthew L. Fidler
