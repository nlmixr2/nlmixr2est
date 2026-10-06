# Calculate Hessian

Unlike \`stats::optimHess\` which assumes the gradient is accurate,
nlmixr2Hess does not make as strong an assumption that the gradient is
accurate but takes more function evaluations to calculate the Hessian.
In addition, this procedures optimizes the forward difference interval
by
[`nlmixr2Gill83`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Gill83.md)

## Usage

``` r
nlmixr2Hess(par, fn, ..., envir = parent.frame())
```

## Arguments

- par:

  Initial values for the parameters to be optimized over.

- fn:

  A function to be minimized (or maximized), with first argument the
  vector of parameters over which minimization is to take place. It
  should return a scalar result.

- ...:

  Extra arguments sent to
  [`nlmixr2Gill83`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Gill83.md)

- envir:

  an environment within which to evaluate the call. This will be most
  useful if `what` is a character string and the arguments are symbols
  or quoted expressions.

## Value

Hessian matrix based on Gill83

## Details

If you have an analytical gradient function, you should use
\`stats::optimHess\`

## See also

[`nlmixr2Gill83`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Gill83.md),
[`optimHess`](https://rdrr.io/r/stats/optim.html)

## Author

Matthew Fidler

## Examples

``` r
 func0 <- function(x){ sum(sin(x))  }
 x <- (0:10)*2*pi/10
 nlmixr2Hess(x, func0)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#>              [,1] [,2]          [,3]          [,4] [,5]        [,6] [,7] [,8]
#>  [1,] 4.15163e-05  NaN  0.000000e+00  0.000000e+00  NaN 0.00000e+00    0  NaN
#>  [2,]         NaN  NaN           NaN           NaN  NaN         NaN  NaN  NaN
#>  [3,] 0.00000e+00  NaN  3.551902e+12 -1.727356e-03  NaN 0.00000e+00    0  NaN
#>  [4,] 0.00000e+00  NaN -1.727356e-03 -2.173232e+12  NaN 0.00000e+00    0  NaN
#>  [5,]         NaN  NaN           NaN           NaN  NaN         NaN  NaN  NaN
#>  [6,] 0.00000e+00  NaN  0.000000e+00  0.000000e+00  NaN 2.15145e-06    0  NaN
#>  [7,] 0.00000e+00  NaN  0.000000e+00  0.000000e+00  NaN 0.00000e+00    0  NaN
#>  [8,]         NaN  NaN           NaN           NaN  NaN         NaN  NaN  NaN
#>  [9,] 0.00000e+00  NaN  0.000000e+00  0.000000e+00  NaN 0.00000e+00    0  NaN
#> [10,] 0.00000e+00  NaN  0.000000e+00  0.000000e+00  NaN 0.00000e+00    0  NaN
#> [11,] 0.00000e+00  NaN  0.000000e+00  0.000000e+00  NaN 0.00000e+00    0  NaN
#>               [,9]    [,10]        [,11]
#>  [1,] 0.0000000000   0.0000 0.000000e+00
#>  [2,]          NaN      NaN          NaN
#>  [3,] 0.0000000000   0.0000 0.000000e+00
#>  [4,] 0.0000000000   0.0000 0.000000e+00
#>  [5,]          NaN      NaN          NaN
#>  [6,] 0.0000000000   0.0000 0.000000e+00
#>  [7,] 0.0000000000   0.0000 0.000000e+00
#>  [8,]          NaN      NaN          NaN
#>  [9,] 0.0001769324   0.0000 0.000000e+00
#> [10,] 0.0000000000 190.1848 0.000000e+00
#> [11,] 0.0000000000   0.0000 3.478511e-06

fr <- function(x) {   ## Rosenbrock Banana function
    x1 <- x[1]
    x2 <- x[2]
    100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { ## Gradient of 'fr'
    x1 <- x[1]
    x2 <- x[2]
    c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
       200 *      (x2 - x1 * x1))
}

h1 <- optimHess(c(1.2,1.2), fr, grr)

h2 <- optimHess(c(1.2,1.2), fr)

## in this case h3 is closer to h1 where the gradient is known

h3 <- nlmixr2Hess(c(1.2,1.2), fr)
```
