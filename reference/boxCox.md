# Cox Box, Yeo Johnson and inverse transformation

Cox Box, Yeo Johnson and inverse transformation

## Usage

``` r
boxCox(x, lambda = 1)

iBoxCox(x, lambda = 1)

yeoJohnson(x, lambda = 1)

iYeoJohnson(x, lambda = 1)
```

## Arguments

- x:

  data to transform

- lambda:

  Cox-box lambda parameter

## Value

Cox-Box Transformed Data

## Author

Matthew L. Fidler

## Examples

``` r
boxCox(1:3,1) ## Normal
#> [1] 0 1 2
iBoxCox(boxCox(1:3,1))
#> [1] 1 2 3

boxCox(1:3,0) ## Log-Normal
#> [1] 0.0000000 0.6931472 1.0986123
iBoxCox(boxCox(1:3,0),0)
#> [1] 1 2 3

boxCox(1:3,0.5) ## lambda=0.5
#> [1] 0.0000000 0.8284271 1.4641016
iBoxCox(boxCox(1:3,0.5),0.5)
#> [1] 1 2 3

yeoJohnson(seq(-3,3),1) ## Normal
#> [1] -3 -2 -1  0  1  2  3
iYeoJohnson(yeoJohnson(seq(-3,3),1))
#> [1] -3 -2 -1  0  1  2  3

yeoJohnson(seq(-3,3),0)
#> [1] -7.5000000 -4.0000000 -1.5000000  0.0000000  0.6931472  1.0986123  1.3862944
iYeoJohnson(yeoJohnson(seq(-3,3),0),0)
#> [1] -3 -2 -1  0  1  2  3
```
