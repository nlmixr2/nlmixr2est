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

A parameter that \`which\` (passed on to
[`nlmixr2Gill83`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Gill83.md))
leaves out of the step search is differenced with the interval the
search would start from, the one \`gillK = 0\` gives every parameter.

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
#>                [,1]          [,2]          [,3]          [,4]          [,5]
#>  [1,] -3.789561e-14  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [2,]  0.000000e+00 -5.877853e-01  0.000000e+00  0.000000e+00 -1.993609e-14
#>  [3,]  0.000000e+00  0.000000e+00 -9.510565e-01 -3.225728e-14  0.000000e+00
#>  [4,]  0.000000e+00  0.000000e+00 -3.225728e-14 -9.510565e-01  0.000000e+00
#>  [5,]  0.000000e+00 -1.993609e-14  0.000000e+00  0.000000e+00 -5.877853e-01
#>  [6,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [7,]  0.000000e+00 -1.994237e-14  0.000000e+00  0.000000e+00  0.000000e+00
#>  [8,]  0.000000e+00  0.000000e+00 -3.226743e-14  0.000000e+00  0.000000e+00
#>  [9,]  0.000000e+00  0.000000e+00  0.000000e+00 -3.226742e-14  0.000000e+00
#> [10,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -1.994236e-14
#> [11,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>                [,6]          [,7]          [,8]          [,9]         [,10]
#>  [1,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>  [2,]  0.000000e+00 -1.994237e-14  0.000000e+00  0.000000e+00  0.000000e+00
#>  [3,]  0.000000e+00  0.000000e+00 -3.226743e-14  0.000000e+00  0.000000e+00
#>  [4,]  0.000000e+00  0.000000e+00  0.000000e+00 -3.226742e-14  0.000000e+00
#>  [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -1.994236e-14
#>  [6,] -5.302314e-14 -1.683612e-13 -2.141588e-13 -2.141587e-13 -1.683612e-13
#>  [7,] -1.683612e-13  5.877853e-01  0.000000e+00  0.000000e+00 -1.993608e-14
#>  [8,] -2.141588e-13  0.000000e+00  9.510565e-01 -3.225726e-14  0.000000e+00
#>  [9,] -2.141587e-13  0.000000e+00 -3.225726e-14  9.510565e-01  0.000000e+00
#> [10,] -1.683612e-13 -1.993608e-14  0.000000e+00  0.000000e+00  5.877853e-01
#> [11,] -1.507584e-14  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00
#>               [,11]
#>  [1,]  0.000000e+00
#>  [2,]  0.000000e+00
#>  [3,]  0.000000e+00
#>  [4,]  0.000000e+00
#>  [5,]  0.000000e+00
#>  [6,] -1.507584e-14
#>  [7,]  0.000000e+00
#>  [8,]  0.000000e+00
#>  [9,]  0.000000e+00
#> [10,]  0.000000e+00
#> [11,]  0.000000e+00

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
