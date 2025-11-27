# Return the objective function

Return the objective function

## Usage

``` r
ofv(x, type, ...)
```

## Arguments

- x:

  object to return objective function value

- type:

  Objective function type value to retrieve or add.

  - focei For most models you can specify "focei" and it will add the
    focei objective function.

  - nlme This switches/chooses the nlme objective function if
    applicable. This objective function cannot be added if it isn't
    present.

  - fo FO objective function value. Cannot be generated

  - foce FOCE object function value. Cannot be generated

  - laplace# This adds/retrieves the Laplace objective function value.
    The `#` represents the number of standard deviations requested when
    expanding the Gaussian Quadrature. This can currently only be used
    with saem fits.

  - gauss#.# This adds/retrieves the Gaussian Quadrature approximation
    of the objective function. The first number is the number of nodes
    to use in the approximation. The second number is the number of
    standard deviations to expand upon.

- ...:

  Other arguments sent to ofv for other methods.

## Value

Objective function value

## Author

Matthew Fidler
