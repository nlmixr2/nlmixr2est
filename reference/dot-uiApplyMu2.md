# This is an internal function for modifying the UI to apply mu2 referencing

mu2 referencing is algebraic mu-referencing by converting to the
transformation to a single value in the original dataset, and moving
that around

## Usage

``` r
.uiApplyMu2(ui, est, data, control)
```

## Arguments

- ui:

  the ui for the model

- est:

  the estimation method

- data:

  the data provided

- control:

  the control object

## Value

Either the original model() block (if changed) or NULL if not changed

## Author

Matthew L. Fidler
