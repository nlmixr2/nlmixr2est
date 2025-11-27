# Augment Prediction for Ipred Model

This function augments the prediction for an individual prediction
(Ipred) model. It retrieves the simulation model from the fit object and
evaluates the model variables.

## Usage

``` r
.augPredIpredModel(fit)
```

## Arguments

- fit:

  The fitted model object from which to retrieve the simulation model.

## Value

The evaluated model variables for the Ipred model.

## Details

The function performs the following steps:

\- Retrieves the simulation model from the provided \`fit\` object using
\`.getSimModel\` with \`hideIpred\` and \`tad\` set to \`FALSE\`.

\- Evaluates the model variables using \`rxModelVars\`.
