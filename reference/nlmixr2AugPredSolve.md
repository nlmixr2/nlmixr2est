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

  specifies the interpolation method for time-varying covariates. When
  solving ODEs it often samples times outside the sampling time
  specified in `events`. When this happens, the time varying covariates
  are interpolated. Currently this can be:

  - `"linear"` interpolation, which interpolates the covariate by
    solving the line between the observed covariates and extrapolating
    the new covariate value.

  - `"locf"` – Last observation carried forward (the default).

  - `"nocb"` – Next Observation Carried Backward. This is the same
    method that NONMEM uses.

  - `"midpoint"` Last observation carried forward to midpoint; Next
    observation carried backward to midpoint.

    For time-varying covariates where a missing value is present, the
    interpolation method will use either "locf" or "nocb" throughout if
    they are the type of covariate interpolation that is selected.

    When using the linear or midpoint interpolation, the lower point in
    the interpolation will use locf to interpolate missing covariates
    and the upper point will use the nocb to interpolate missing
    covariates.

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
