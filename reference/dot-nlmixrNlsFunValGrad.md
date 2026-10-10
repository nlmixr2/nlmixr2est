# The residuals and their Jacobian for nls

The right-hand side of the
[`stats::nls()`](https://rdrr.io/r/stats/nls.html) formula
(`ui$nlsFormula`).

## Usage

``` r
.nlmixrNlsFunValGrad(DV, ...)
```

## Arguments

- DV:

  dependent variable

- ...:

  The estimated parameters (scaled)

## Value

The residuals of the loaded nls problem, with their Jacobian as the
`"gradient"` attribute

## Details

This is an internal function and should not be called directly.

## Author

Matthew L. Fidler
