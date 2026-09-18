# The cache key of a covariance method's options

[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md)
records the result of this generic for every covariance it computes, and
reuses a cached covariance only when the key it asks for is
[`identical()`](https://rdrr.io/r/base/identical.html) to the recorded
one. The default key is the control object itself, as a plain list.

## Usage

``` r
setCovOptions(control, fit, ...)

# Default S3 method
setCovOptions(control, fit, ...)

# S3 method for class 'rsControl'
setCovOptions(control, fit, ...)
```

## Arguments

- control:

  covariance control object (for example
  [`rsControl()`](https://nlmixr2.github.io/nlmixr2est/reference/rsControl.md))

- fit:

  nlmixr2 fit, or its environment when the key of a built-in method's
  defaults is taken

- ...:

  ignored

## Value

named list, compared with
[`identical()`](https://rdrr.io/r/base/identical.html)

## Details

A package whose covariance depends on more than its control – a
covariance seeded from another covariance on the fit, say – adds a
method for its control class that puts that state in the key, so a
change in it recomputes the covariance instead of reinstalling a stale
one. Options that do not change the result (parallel workers, say) can
be left out.

## See also

[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md)

## Author

Matt Fidler

## Examples

``` r
setCovOptions(saControl(), NULL)
#> $nBurn
#> [1] 100
#> 
#> $nEm
#> [1] 100
#> 
#> $nSaCov
#> [1] 500
#> 
#> $seed
#> [1] 99
#> 
```
