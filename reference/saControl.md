# Options for the SAEM stochastic-approximation covariance in setCov()

Used by `setCov(fit, "sa")`, which runs a short SAEM at the fit's
estimates before the covariance phase.

## Usage

``` r
saControl(nBurn = 100L, nEm = 100L, nSaCov = 500L, seed = 99L)
```

## Arguments

- nBurn, nEm:

  warm-up iterations that equilibrate the MCMC chains before the
  covariance phase

- nSaCov:

  iterations in the covariance phase; more gives a less noisy covariance

- seed:

  random seed

## Value

`saControl` object

## See also

[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Author

Matt Fidler

## Examples

``` r
saControl(nSaCov = 1000)
#> $nBurn
#> [1] 100
#> 
#> $nEm
#> [1] 100
#> 
#> $nSaCov
#> [1] 1000
#> 
#> $seed
#> [1] 99
#> 
#> attr(,"class")
#> [1] "saControl"
```
