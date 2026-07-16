# Get the ODE states of a model (rxode2 v3/v4 compatible)

Calls
[`rxode2::rxState()`](https://nlmixr2.github.io/rxode2/reference/rxState.html)
(or
[`rxode2::rxStateOde()`](https://nlmixr2.github.io/rxode2/reference/rxStateOde.html)
with rxode2 version 4) on the input.

## Usage

``` r
rxode2stateOde(inp)
```

## Arguments

- inp:

  rxode2 model (or symengine environment) to query

## Value

character vector of ODE state names

## Author

Matthew L. Fidler
