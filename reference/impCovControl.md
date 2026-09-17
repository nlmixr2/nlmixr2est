# Options for the importance-sampling covariance in setCov()

Used by `setCov(fit, "imp")`, which runs frozen importance-sampling EM
iterations at the fit's estimates.

## Usage

``` r
impCovControl(nIter = 1L, isample = 300L, impSeed = 42L)
```

## Arguments

- nIter:

  frozen EM iterations (`0` is an E-step-only evaluation)

- isample:

  Number of importance samples drawn per subject per iteration (NONMEM
  ISAMPLE). Either a single count used for every subject, or a vector of
  length \`nsub\` giving a count \*\*per subject\*\*.

  Per-subject counts are the NM7 Technical Guide's own remedy for poor
  coverage (its derivation is Gaussian throughout and never mentions a t
  proposal): a subject whose weights are badly behaved can be given more
  samples without charging every other subject for them. Note this
  treats the symptom rather than the cause – more draws from a proposal
  whose tails are too light still gives weights with infinite variance,
  which \`fit\$env\$impPsisK\` will show. See \`df\` for the shape-based
  remedy.

- impSeed:

  Base seed for the per-subject thread-safe (threefry) RNG streams;
  results are reproducible and independent of the thread count.

## Value

`impCovControl` object

## See also

[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md),
[`impmapControl()`](https://nlmixr2.github.io/nlmixr2est/reference/impmapControl.md)

## Author

Matt Fidler

## Examples

``` r
impCovControl(isample = 1000)
#> $nIter
#> [1] 1
#> 
#> $isample
#> [1] 1000
#> 
#> $impSeed
#> [1] 42
#> 
#> attr(,"class")
#> [1] "impCovControl"
```
