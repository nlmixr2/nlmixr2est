# Options for the SAEM stochastic-approximation covariance in setCov()

Used by `setCov(fit, "sa")`, which runs the SAEM covariance phase at the
fit's estimates, after warm-up iterations unless it continues a SAEM
fit's own chains. Every population parameter is held at the fit's
estimates throughout (mixture proportions excepted), so the covariance
is the one at those estimates.

## Usage

``` r
saControl(
  nBurn = 100L,
  nEm = 100L,
  nSaCov = 500L,
  seed = 99L,
  warmStart = TRUE
)
```

## Arguments

- nBurn, nEm:

  warm-up iterations that equilibrate the MCMC chains before the
  covariance phase, for a fit without SAEM chain history (or with
  `warmStart = FALSE`)

- nSaCov:

  iterations in the covariance phase; more gives a less noisy covariance

- seed:

  random seed

- warmStart:

  for a SAEM fit, continue the fit's own MCMC chains: the covariance
  phase starts from the fit's last iteration with no new warm-up
  iterations. This changes the Monte Carlo path, so the covariance
  differs from a cold start by Monte Carlo noise. Other fits, and SAEM
  fits that kept no chain history, always run `nBurn`/`nEm`.

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
#> $warmStart
#> [1] TRUE
#> 
#> attr(,"class")
#> [1] "saControl"
```
