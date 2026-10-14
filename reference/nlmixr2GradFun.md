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

  Max steps to determine the optimal forward/central difference step
  size per parameter (Gill 1983). \`0\` = no optimal step size
  determined.

- gillStep:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration the new step size = (prior step size)\*gillStep

- gillFtol:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates.

- useColor:

  Logical (or \`NULL\`) emit ANSI bold/color escapes in the iteration
  print. \`NULL\` (default) defers to \[crayon::has_color()\].

- printNcol:

  Integer (or \`NULL\`) parameter columns per row before wrapping.
  \`NULL\` (default) uses \`floor((getOption("width") - 23) / 12)\`.

- print:

  Either a scalar print-frequency (\`0\` = suppress, \`1\` (default) =
  every evaluation, \`N\` = every Nth), OR a pre-built
  \[iterPrintControl()\] object. Equivalent to \`iterPrintControl(every
  = print, ncol = printNcol, useColor = useColor)\`.

## Value

A list with \`eval\`, \`grad\`, \`hist\` and \`unscaled\` functions.
This is an internal module used with dynmodel

## Examples

``` r

func0 <- function(x){ sum(sin(x))  }

## This will printout every iteration or when print=X
gf <- nlmixr2GradFun(func0)

## x
x <- (0:10)*2*pi/10;
gf$eval(x)
#> [1] -1.224145e-16
gf$grad(x)
#>  [1]  1.0000000  0.8089234  0.3088979 -0.3091360 -0.8091106 -0.9999993
#>  [7] -0.8089234 -0.3088979  0.3091360  0.8091106  0.9999979

## x2
x2 <- x+0.1
gf$eval(x2)
#> [1] 0.09983342
gf$grad(x2)
#>  [1]  0.9949773  0.7461822  0.2123982 -0.4025390 -0.8637384 -0.9949001
#>  [7] -0.7461865 -0.2124014  0.4025369  0.8637373  0.9948221

## Gives the parameter history as a data frame
gf$hist()
#>   iter               type          objf        t1        t2        t3        t4
#> 1    1           Unscaled -1.224145e-16 0.0000000 0.6283185 1.2566371  1.884956
#> 2    2           Unscaled  9.983342e-02 0.1000000 0.7283185 1.3566371  1.984956
#> 3    1    Gill83 Gradient            NA 1.0000000 0.8089234 0.3088979 -0.309136
#> 4    2 Forward Difference            NA 0.9949773 0.7461822 0.2123982 -0.402539
#>           t5         t6         t7         t8        t9       t10       t11
#> 1  2.5132741  3.1415927  3.7699112  4.3982297 5.0265482 5.6548668 6.2831853
#> 2  2.6132741  3.2415927  3.8699112  4.4982297 5.1265482 5.7548668 6.3831853
#> 3 -0.8091106 -0.9999993 -0.8089234 -0.3088979 0.3091360 0.8091106 0.9999979
#> 4 -0.8637384 -0.9949001 -0.7461865 -0.2124014 0.4025369 0.8637373 0.9948221
```
