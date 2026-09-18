# Install an already-computed covariance result on a fit

`setCov(fit) <- value` installs the covariance carried by `value` as the
fit's covariance, the way `setCov(fit, method)` installs one it
computes: the matrix is checked for positive definiteness, the standard
errors are refreshed, and the prior covariance is kept in `fit$covList`.
The options that produced `value` are recorded too, so a later
`setCov(fit, method, control = ...)` reuses it only when the options it
asks for are the same.

## Usage

``` r
setCov(fit, method = NULL, ...) <- value
```

## Arguments

- fit:

  nlmixr2 fit

- method:

  covariance-method name to install `value` as; `NULL` lets
  [`setCovValue()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovValue.md)
  choose

- ...:

  passed to
  [`setCovValue()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovValue.md)

- value:

  a covariance matrix named like `fit$cov`, or a result with a
  [`setCovValue()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovValue.md)
  method

## Value

the fit, with its covariance updated

## Details

`value` is dispatched on through
[`setCovValue()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovValue.md):
a covariance matrix is installed as `method` (default `"user"`) with no
options, and other packages add methods for their own results (for
example a SIR run), returning the matrix, the method name and the
options.

## See also

[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md),
[`setCovOptions()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovOptions.md)

## Author

Matt Fidler

## Examples

``` r
if (FALSE) { # \dontrun{
setCov(fit) <- fit$cov * 2
setCov(fit, "doubled") <- fit$cov * 2
} # }
```
