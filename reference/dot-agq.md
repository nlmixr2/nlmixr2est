# Get the adaptive Gauss-Hermite quadrature points and weights

Get the adaptive Gauss-Hermite quadrature points and weights

## Usage

``` r
.agq(neta = 2, nAGQ = 3)
```

## Arguments

- neta:

  number of eta parameters in the model

- nAGQ:

  number of adaptive quadrature points to use

## Value

A list with the following elements:

- x:

  A matrix of quadrature points, one row per point.

- w:

  A matrix of quadrature weights, one row per point.

- n:

  The number of quadrature points.

- neta:

  The number of eta parameters.

- nAQD:

  The number of adaptive quadrature points.

- first:

  A logical indicating if the first point is zero.

## Author

Matthew L. Fidler

## Examples

``` r
.agq(neta=2, nAGQ=3)
#> $x
#>            Var1      Var2
#>  [1,]  0.000000  0.000000
#>  [2,] -1.224745  0.000000
#>  [3,]  1.224745  0.000000
#>  [4,]  0.000000 -1.224745
#>  [5,] -1.224745 -1.224745
#>  [6,]  1.224745 -1.224745
#>  [7,]  0.000000  1.224745
#>  [8,] -1.224745  1.224745
#>  [9,]  1.224745  1.224745
#> 
#> $w
#>            Var1      Var2
#>  [1,] 0.6666667 0.6666667
#>  [2,] 0.1666667 0.6666667
#>  [3,] 0.1666667 0.6666667
#>  [4,] 0.6666667 0.1666667
#>  [5,] 0.1666667 0.1666667
#>  [6,] 0.1666667 0.1666667
#>  [7,] 0.6666667 0.1666667
#>  [8,] 0.1666667 0.1666667
#>  [9,] 0.1666667 0.1666667
#> 
#> $n
#> [1] 9
#> 
#> $neta
#> [1] 2
#> 
#> $nAGQ
#> [1] 3
#> 
#> $first
#> [1] TRUE
#> 
```
