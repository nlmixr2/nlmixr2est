# This downgrades the UI for any of the zero etas in the model

This downgrades the UI for any of the zero etas in the model

## Usage

``` r
.downgradeEtas(ui, zeroEtas = character(0))
```

## Arguments

- ui:

  rxode2 User interface function

- zeroEtas:

  The names of the zero etas in the model

## Value

New rxode2 ui with the zero etas removed

## Author

Matthew L. Fidler
