# Create a gradient function based on gill numerical differences

Create a gradient function based on gill numerical differences

## Usage

``` r
nlmixr2Eval_(theta, md5)

nlmixr2Unscaled_(theta, md5)

nlmixr2Grad_(theta, md5)

nlmixr2ParHist_(md5)

nlmixr2GradFun(
  what,
  envir = parent.frame(),
  which,
  thetaNames,
  gillRtol = sqrt(.Machine$double.eps),
  gillK = 10L,
  gillStep = 2,
  gillFtol = 0,
  useColor = crayon::has_color(),
  printNcol = floor((getOption("width") - 23)/12),
  print = 1
)
```

## Arguments

- theta:

  for the internal functions theta is the parameter values

- md5:

  the md5 identifier for the internal gradient function information.

- what:

  either a function or a non-empty character string naming the function
  to be called.

- envir:

  an environment within which to evaluate the call. This will be most
  useful if `what` is a character string and the arguments are symbols
  or quoted expressions.

- which:

  Which parameters to calculate the forward difference and optimal
  forward difference interval

- thetaNames:

  Names for the theta parameters

- gillRtol:

  The relative tolerance used for Gill 1983 determination of optimal
  step size.

- gillK:

  The total number of possible steps to determine the optimal
  forward/central difference step size per parameter (by the Gill 1983
  method). If 0, no optimal step size is determined. Otherwise this is
  the optimal step size determined.

- gillStep:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration the new step size = (prior step size)\*gillStep

- gillFtol:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates.

- useColor:

  Boolean indicating if focei can use ASCII color codes

- printNcol:

  Number of columns to printout before wrapping parameter
  estimates/gradient

- print:

  Integer representing when the outer step is printed. When this is 0 or
  do not print the iterations. 1 is print every function evaluation
  (default), 5 is print every 5 evaluations.

## Value

A list with \`eval\`, \`grad\`, \`hist\` and \`unscaled\` functions.
This is an internal module used with dynmodel

## Examples

``` r
func0 <- function(x){ sum(sin(x))  }

## This will printout every interation or when print=X
gf <- nlmixr2GradFun(func0)

## x
x <- (0:10)*2*pi/10;
gf$eval(x)
#> [1] -2.291621e-08
gf$grad(x)
#>  [1]  1.0000000  0.8090168  0.3090166 -0.3090170 -0.8090170 -1.0000000
#>  [7] -0.8090170 -0.3090170  0.3090170  0.8090170  1.0000000

## x2
x2 <- x+0.1
gf$eval(x2)
#> [1] 0.1001453
gf$grad(x2)
#>  [1]  0.9950042  0.7462946  0.2125259 -0.4024205 -0.8636559 -0.9950354
#>  [7] -0.7462947 -0.2125260  0.4024204  0.8636559  0.9949954

## Gives the parameter history as a data frame
gf$hist()
#>   iter               type          objf        t1        t2        t3
#> 1    1           Unscaled -2.291621e-08 0.0000000 0.6283185 1.2566371
#> 2    2           Unscaled  1.001453e-01 0.1000000 0.7283185 1.3566371
#> 3    1    Gill83 Gradient            NA 1.0000000 0.8090168 0.3090166
#> 4    2 Forward Difference            NA 0.9950042 0.7462946 0.2125259
#>           t4         t5         t6         t7        t8        t9       t10
#> 1  1.8849556  2.5132741  3.1415926  3.7699112  4.398230 5.0265482 5.6548667
#> 2  1.9849556  2.6132741  3.2412792  3.8699112  4.498230 5.1265482 5.7548667
#> 3 -0.3090170 -0.8090170 -1.0000000 -0.8090170 -0.309017 0.3090170 0.8090170
#> 4 -0.4024205 -0.8636559 -0.9950354 -0.7462947 -0.212526 0.4024204 0.8636559
#>         t11
#> 1 6.2831853
#> 2 6.3831853
#> 3 1.0000000
#> 4 0.9949954
```
