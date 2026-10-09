# Update the nlmixr2 object with new fit information

Update the nlmixr2 object with new fit information

## Usage

``` r
nlmixrUpdateObject(fit, objName, envir, origFitEnv = NULL)
```

## Arguments

- fit:

  nlmixr2 fit object to update in the environment

- objName:

  Name of the object

- envir:

  Environment to search

- origFitEnv:

  Original fit\$env to compare, otherwise simply use fit\$env

## Value

`TRUE` (invisibly) when the binding was updated, `FALSE` otherwise (e.g.
`objName` is not a single name, as for `addCwres(fits[[1]])`)

## Author

Matthew L. Fidler
