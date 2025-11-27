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
#>               [,1]        [,2]        [,3]        [,4]        [,5]       [,6]
#>  [1,] 143438578824           0           0           0           0          0
#>  [2,]            0 54098751870           0           0           0          0
#>  [3,]            0           0 28167126506           0           0          0
#>  [4,]            0           0           0 17234064451           0          0
#>  [5,]            0           0           0           0 11620957271          0
#>  [6,]            0           0           0           0           0 8362405349
#>  [7,]            0           0           0           0           0          0
#>  [8,]            0           0           0           0           0          0
#>  [9,]            0           0           0           0           0          0
#> [10,]            0           0           0           0           0          0
#> [11,]            0           0           0           0           0          0
#>               [,7]         [,8]         [,9]        [,10]    [,11]
#>  [1,]            0            0            0            0     0.00
#>  [2,]            0            0            0            0     0.00
#>  [3,]            0            0            0            0     0.00
#>  [4,]            0            0            0            0     0.00
#>  [5,]            0            0            0            0     0.00
#>  [6,]            0            0            0            0     0.00
#>  [7,] 337177531807            0            0            0     0.00
#>  [8,]            0 546113920621            0            0     0.00
#>  [9,]            0            0 545713835656            0     0.00
#> [10,]            0            0            0 337491323869     0.00
#> [11,]            0            0            0            0 26296.35

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
