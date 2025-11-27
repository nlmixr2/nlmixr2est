# Add objective function data frame to the current objective function

Add objective function data frame to the current objective function

## Usage

``` r
nlmixrAddObjectiveFunctionDataFrame(fit, objDf, type, etaObf = NULL)
```

## Arguments

- fit:

  nlmixr fit object

- objDf:

  nlmixr objective function data frame which has column names "OBJF",
  "AIC", "BIC", "Log-likelihood" and "Condition#(Cov)" "Condition#(Cor)"

- type:

  Objective Function Type

- etaObf:

  Eta objective function table to add (with focei) to give focei
  objective function

## Value

Nothing, called for side effects

## Author

Matthew L. Fidler
