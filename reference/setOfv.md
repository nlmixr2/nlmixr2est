# Set/get Objective function type for a nlmixr2 object

Set/get Objective function type for a nlmixr2 object

## Usage

``` r
setOfv(x, type)

getOfvType(x)
```

## Arguments

- x:

  nlmixr2 fit object

- type:

  Type of objective function to use for AIC, BIC, and \$objective.
  `"imp"` and `"impmap"` add an importance-sampling objective from an
  E-step-only run (`nIter=0`) at the fit's estimates.

## Value

Nothing

## Author

Matthew L. Fidler
