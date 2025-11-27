# This function process the data for use in focei

The \$origData is the data that is fed into the focei before
modification The \$dataSav is the data saved for focei

## Usage

``` r
.foceiPreProcessData(data, env, ui, rxControl = NULL)
```

## Arguments

- data:

  Input dataset

- env:

  focei environment where focei family is run

- ui:

  rxode2 ui

- rxControl:

  is the rxode2 control that is used to translate to the modeling
  dataset

## Value

Nothing, called for side effects

## Author

Matthew L. Fidler
