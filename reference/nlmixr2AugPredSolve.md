# Augmented Prediction for nlmixr2 fit

Augmented Prediction for nlmixr2 fit

## Usage

``` r
nlmixr2AugPredSolve(
  fit,
  covsInterpolation = c("locf", "nocb", "linear", "midpoint"),
  minimum = NULL,
  maximum = NULL,
  length.out = 51L,
  ...
)

# S3 method for class 'nlmixr2FitData'
augPred(
  object,
  primary = NULL,
  minimum = NULL,
  maximum = NULL,
  length.out = 51,
  ...
)
```

## Arguments

- fit:

  Nlmixr2 fit object

- covsInterpolation:

  covariate interpolation; by default the one the fit was estimated with

- minimum:

  an optional lower limit for the primary covariate. Defaults to
  `min(primary)`.

- maximum:

  an optional upper limit for the primary covariate. Defaults to
  `max(primary)`.

- length.out:

  an optional integer with the number of primary covariate values at
  which to evaluate the predictions. Defaults to 51.

- ...:

  some methods for the generic may require additional arguments.

- object:

  a fitted model object from which predictions can be extracted, using a
  `predict` method.

- primary:

  an optional one-sided formula specifying the primary covariate to be
  used to generate the augmented predictions. By default, if a covariate
  can be extracted from the data used to generate `object` (using
  `getCovariate`), it will be used as `primary`.

## Value

Stacked data.frame with observations, individual/population predictions.

## Author

Matthew L. Fidler
