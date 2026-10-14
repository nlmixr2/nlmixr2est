# Get the optimal forward difference interval by Gill83 method

Get the optimal forward difference interval by Gill83 method

## Usage

``` r
nlmixr2Gill83(
  what,
  args,
  envir = parent.frame(),
  which,
  gillRtol = sqrt(.Machine$double.eps),
  gillK = 10L,
  gillStep = 2,
  gillFtol = 0
)
```

## Arguments

- what:

  either a function or a non-empty character string naming the function
  to be called.

- args:

  a *list* of arguments to the function call. The `names` attribute of
  `args` gives the argument names.

- envir:

  an environment within which to evaluate the call. This will be most
  useful if `what` is a character string and the arguments are symbols
  or quoted expressions.

- which:

  Which parameters to calculate the forward difference and optimal
  forward difference interval

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

## Value

A data frame with the following columns:

\- info Gradient evaluation/forward difference information

\- hf Forward difference final estimate

\- df Derivative estimate

\- df2 2nd Derivative Estimate

\- err Error of the final estimate derivative

\- aEps Absolute difference for forward numerical differences

\- rEps Relative Difference for backward numerical differences

\- aEpsC Absolute difference for central numerical differences

\- rEpsC Relative difference for central numerical differences

The `info` returns one of the following:

\- "Not Assessed" Gradient wasn't assessed: the parameter is left out by
`which`, or `gillK = 0`, which determines no interval and reports the
one the search starts from as `hf`

\- "Good Success" in Estimating optimal forward difference interval

\- "High Grad Error" Large error; Derivative estimate error `fTol` or
more of the derivative

\- "Constant Grad" Function constant or nearly constant for this
parameter

\- "Odd/Linear Grad" Function odd or nearly linear, df = K, df2 ~ 0

\- "Grad changes quickly" df2 increases rapidly as h decreases

## Author

Matthew Fidler

## Examples

``` r

## These are taken from the numDeriv's grad examples to show how
## simple gradients are assessed with nlmixr2Gill83

nlmixr2Gill83(sin, pi)
#> Gill83 Derivative/Forward Difference
#>   (rtol=1.49011611938477e-08; K=10, step=2, ftol=0)
#> 
#>   info          hf      hphi         df df2          err         aEps
#> 1 Good 0.002022262 0.0161781 -0.9999993   0 1.473712e-05 0.0004882812
#>           rEps        aEpsC        rEpsC            f
#> 1 0.0004882812 0.0004882812 0.0004882812 1.224647e-16

nlmixr2Gill83(sin, (0:10)*2*pi/10)
#> Gill83 Derivative/Forward Difference
#>   (rtol=1.49011611938477e-08; K=10, step=2, ftol=0)
#> 
#>                    info           hf         hphi         df           df2
#> 1       Odd/Linear Grad 0.0004882812 0.0002441406  1.0000000             0
#> 2  Grad changes quickly 0.0007950774 0.0003975387  0.8090170 1949978254364
#> 3  Grad changes quickly 0.0011018736 0.0005509368  0.3090170 1642754654097
#> 4  Grad changes quickly 0.0014086697 0.0007043349 -0.3090170 1005119906032
#> 5  Grad changes quickly 0.0017154659 0.0008577329 -0.8090170  418874984766
#> 6       Odd/Linear Grad 0.0020222620 0.0010111310 -0.9999993             0
#> 7  Grad changes quickly 0.0023290582 0.0011645291 -0.8090170 -227241639381
#> 8  Grad changes quickly 0.0026358544 0.0013179272 -0.3090170 -287073706156
#> 9  Grad changes quickly 0.0029426505 0.0014713253  0.3090170 -230334431768
#> 10 Grad changes quickly 0.0032494467 0.0016247233  0.8090170 -116742709124
#> 11                 Good 0.0035562428 0.0142249713  0.9999979             0
#>             err         aEps         rEps        aEpsC        rEpsC
#> 1  6.103516e-05 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 2  7.751918e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 3  9.050540e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 4  7.079410e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 5  3.592829e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 6  1.473712e-05 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 7  2.646295e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 8  3.783422e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 9  3.388969e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 10 1.896746e+08 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#> 11 8.380283e-06 0.0004882812 0.0004882812 0.0004882812 0.0004882812
#>                f
#> 1  -2.449294e-16
#> 2  -2.449294e-16
#> 3  -2.449294e-16
#> 4  -2.449294e-16
#> 5  -2.449294e-16
#> 6  -2.449294e-16
#> 7  -2.449294e-16
#> 8  -2.449294e-16
#> 9  -2.449294e-16
#> 10 -2.449294e-16
#> 11 -2.449294e-16

func0 <- function(x){ sum(sin(x))  }
nlmixr2Gill83(func0 , (0:10)*2*pi/10)
#> Gill83 Derivative/Forward Difference
#>   (rtol=1.49011611938477e-08; K=10, step=2, ftol=0)
#> 
#>    info           hf        hphi         df        df2          err
#> 1  Good 0.0004882812 0.125000000  1.0000000  0.0000000 6.103516e-05
#> 2  Good 0.0003184423 0.001590155  0.8089234 -0.5877851 1.871756e-04
#> 3  Good 0.0002503438 0.001101874  0.3088979 -0.9510564 2.380911e-04
#> 4  Good 0.0002503438 0.001408670 -0.3091360 -0.9510564 2.380911e-04
#> 5  Good 0.0003184423 0.001715466 -0.8091106 -0.5877851 1.871756e-04
#> 6  Good 0.0020222620 0.016178096 -0.9999993  0.0000000 1.473712e-05
#> 7  Good 0.0003184423 0.002329058 -0.8089234  0.5877850 1.871756e-04
#> 8  Good 0.0002503439 0.002635854 -0.3088979  0.9510560 2.380911e-04
#> 9  Good 0.0002503439 0.002942651  0.3091360  0.9510558 2.380910e-04
#> 10 Good 0.0003184424 0.003249447  0.8091106  0.5877847 1.871756e-04
#> 11 Good 0.0035562428 0.028449943  0.9999979  0.0000000 8.380283e-06
#>            aEps         rEps        aEpsC        rEpsC             f
#> 1  4.882812e-04 4.882812e-04 4.882812e-04 4.882812e-04 -1.224145e-16
#> 2  1.955651e-04 1.955651e-04 1.955651e-04 1.955651e-04 -1.224145e-16
#> 3  1.109367e-04 1.109367e-04 1.109367e-04 1.109367e-04 -1.224145e-16
#> 4  8.677563e-05 8.677563e-05 8.677563e-05 8.677563e-05 -1.224145e-16
#> 5  9.063976e-05 9.063976e-05 9.063976e-05 9.063976e-05 -1.224145e-16
#> 6  4.882812e-04 4.882812e-04 4.882812e-04 4.882812e-04 -1.224145e-16
#> 7  6.676064e-05 6.676064e-05 6.676064e-05 6.676064e-05 -1.224145e-16
#> 8  4.637518e-05 4.637518e-05 4.637518e-05 4.637518e-05 -1.224145e-16
#> 9  4.154018e-05 4.154018e-05 4.154018e-05 4.154018e-05 -1.224145e-16
#> 10 4.785106e-05 4.785106e-05 4.785106e-05 4.785106e-05 -1.224145e-16
#> 11 4.882812e-04 4.882812e-04 4.882812e-04 4.882812e-04 -1.224145e-16

func1 <- function(x){ sin(10*x) - exp(-x) }
curve(func1,from=0,to=5)


x <- 2.04
numd1 <- nlmixr2Gill83(func1, x)
exact <- 10*cos(10*x) + exp(-x)
c(numd1$df, exact, (numd1$df - exact)/exact)
#> [1]  0.332315771  0.333537144 -0.003661882

x <- c(1:10)
numd1 <- nlmixr2Gill83(func1, x)
exact <- 10*cos(10*x) + exp(-x)
cbind(numd1=numd1$df, exact, err=(numd1$df - exact)/exact)
#>           numd1     exact           err
#>  [1,] -8.022836 -8.022836 -1.166840e-10
#>  [2,]  4.216156  4.216156  6.611818e-11
#>  [3,]  1.592302  1.592302 -2.766382e-10
#>  [4,] -6.651065 -6.651065 -9.572443e-11
#>  [5,]  9.656398  9.656398 -4.522570e-12
#>  [6,] -9.521651 -9.521651 -2.119119e-10
#>  [7,]  6.334104  6.334104 -3.813201e-10
#>  [8,] -1.103537 -1.103537  7.983780e-11
#>  [9,] -4.480613 -4.480613 -2.876500e-10
#> [10,]  8.624103  8.623234  1.007135e-04

sc2.f <- function(x){
  n <- length(x)
   sum((1:n) * (exp(x) - x)) / n
}

sc2.g <- function(x){
  n <- length(x)
  (1:n) * (exp(x) - 1) / n
}

x0 <- rnorm(100)
exact <- sc2.g(x0)

g <- nlmixr2Gill83(sc2.f, x0)

max(abs(exact - g$df)/(1 + abs(exact)))
#> [1] 0.001120989
```
