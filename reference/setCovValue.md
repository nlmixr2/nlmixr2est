# The covariance, method name and options a result installs

The generic behind `setCov(fit) <- value`. A method returns
`list(cov =, method =, options =)`: the covariance named like `fit$cov`,
the covariance-method name to install it as, and the options recorded
for it (the same key
[`setCovOptions()`](https://nlmixr2.github.io/nlmixr2est/reference/setCovOptions.md)
gives the control that would compute it, so
[`setCov()`](https://nlmixr2.github.io/nlmixr2est/reference/setCov.md)
can reuse it). It may also return `extra`, a named list of objects
stored in the fit environment (so `fit$<name>` returns them) once the
covariance is installed – the full result the covariance came from, say.

## Usage

``` r
setCovValue(value, fit, method = NULL, ...)

# Default S3 method
setCovValue(value, fit, method = NULL, ...)

# S3 method for class 'matrix'
setCovValue(value, fit, method = NULL, ...)
```

## Arguments

- value:

  the result to install

- fit:

  nlmixr2 fit

- method:

  requested covariance-method name, or `NULL`

- ...:

  ignored by the built-in methods

## Value

`list(cov, method, options)`, optionally with `extra`

## See also

`setCov<-`

## Author

Matt Fidler
