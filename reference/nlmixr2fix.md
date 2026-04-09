# Try to fix a nlmixr2 fit

Currently this re-evaluates the function in the current version of
rxode2.

## Usage

``` r
nlmixr2fix(fit)
```

## Arguments

- fit:

  nlmixr2 fit object from a different version of nlmixr2.

## Value

A nlmixr2 fit that has been (possibly) adjusted to work with the current
version of nlmixr2.

## Author

Matthew L. Fidler

## Examples

``` r
if (FALSE) { # \dontrun{


  # This is a nlmixr2 v3 fit and requires the qs package to read in
  fit <- system.file("testfit_nlmixr3.rds", package = "nlmixr2est")
  fit <- readRDS(fit)

  # While it prints well, it can't be used in all functions because
  # Language features (like +var()) are not supported in the v3 version

  print(fit)

  try(rxSolve(fit)) # should error, but with try it will just display the error

  # This function attempts to fix it by regenerating the rxode2 model with the
  # new features

  # This function also prints out the information on how this fit was created

  fit <- try(nlmixr2fix(fit))

  # Now solving and other functions work
  if (!inherits(fit, "try-error")) {
    rxSolve(fit)
  }

} # }
```
