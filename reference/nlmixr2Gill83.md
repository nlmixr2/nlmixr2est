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
#>              info           hf         hphi df df2          err         aEps
#> 1 Odd/Linear Grad 2.237911e-11 1.118956e-11 -1   0 1.630865e-13 5.403504e-12
#>           rEps        aEpsC        rEpsC            f
#> 1 5.403504e-12 5.403504e-12 5.403504e-12 1.224647e-16

nlmixr2Gill83(sin, (0:10)*2*pi/10)
#> Warning: NaNs produced
#> Gill83 Derivative/Forward Difference
#>   (rtol=1.49011611938477e-08; K=10, step=2, ftol=0)
#> 
#>                    info           hf         hphi         df           df2
#> 1  Grad changes quickly 7.641709e-12 3.820854e-12  1.0000000  8.796093e+12
#> 2  Grad changes quickly 1.244314e-11 6.221568e-12  0.8040129  7.961390e+27
#> 3  Grad changes quickly 1.724456e-11 8.622281e-12  0.3098531  6.707055e+27
#> 4  Grad changes quickly 2.204599e-11 1.102299e-11 -0.3094082  4.103713e+27
#> 5  Grad changes quickly 2.684742e-11 1.342371e-11 -0.8087998  1.710187e+27
#> 6  Grad changes quickly 3.164884e-11 1.582442e-11 -1.0057971  7.692125e+11
#> 7  Grad changes quickly 3.645027e-11 1.822514e-11 -0.8078099 -9.277844e+26
#> 8  Grad changes quickly 4.125170e-11 2.062585e-11 -0.3059084 -1.172067e+27
#> 9  Grad changes quickly 4.605313e-11 2.302656e-11  0.3110439 -9.404117e+26
#> 10 Grad changes quickly 5.085455e-11 2.542728e-11  0.8103793 -4.766383e+26
#> 11      High Grad Error          Inf 5.699172e-08        NaN  0.000000e+00
#>             err         aEps         rEps        aEpsC        rEpsC
#> 1  3.360859e+01 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 2  4.953233e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 3  5.783011e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 4  4.523521e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 5  2.295705e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 6  1.217234e+01 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 7  1.690900e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 8  2.417488e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 9  2.165445e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 10 1.211961e+16 7.641709e-12 7.641709e-12 7.641709e-12 7.641709e-12
#> 11          NaN          Inf          Inf          Inf          Inf
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
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Gill83 Derivative/Forward Difference
#>   (rtol=1.49011611938477e-08; K=10, step=2, ftol=0)
#> 
#>               info           hf         hphi         df           df2
#> 1  Odd/Linear Grad 5.402397e-12 2.701199e-12  1.0000000  0.000000e+00
#> 2  High Grad Error          Inf 8.796823e-12        NaN  0.000000e+00
#> 3             Good 6.250714e-15 2.438250e-11  0.3197077  1.867471e+05
#> 4  High Grad Error 7.991108e-15 3.117135e-11 -0.3195443 -1.142612e+05
#> 5  High Grad Error          Inf 4.745025e-12        NaN  0.000000e+00
#> 6  Odd/Linear Grad 2.237453e-11 1.118726e-11 -1.0000000  0.000000e+00
#> 7  High Grad Error 1.321229e-14 5.153791e-11 -0.8150866  4.179811e+04
#> 8  High Grad Error          Inf 1.458169e-11        NaN  0.000000e+00
#> 9             Good 2.670893e-13 1.041850e-09  0.3092621  1.022822e+02
#> 10            Good 8.542248e-10 1.797612e-11  0.8090174  9.999278e-06
#> 11 Odd/Linear Grad 3.934666e-11 1.967333e-11  1.0000000  0.000000e+00
#>             err         aEps         rEps        aEpsC        rEpsC
#> 1  6.752996e-13 5.402397e-12 5.402397e-12 5.402397e-12 5.402397e-12
#> 2           NaN          Inf          Inf          Inf          Inf
#> 3  1.167302e-09 2.769924e-15 2.769924e-15 2.769924e-15 2.769924e-15
#> 4  9.130740e-10 2.769924e-15 2.769924e-15 2.769924e-15 2.769924e-15
#> 5           NaN          Inf          Inf          Inf          Inf
#> 6  1.630531e-13 5.402397e-12 5.402397e-12 5.402397e-12 5.402397e-12
#> 7  5.522489e-10 2.769924e-15 2.769924e-15 2.769924e-15 2.769924e-15
#> 8           NaN          Inf          Inf          Inf          Inf
#> 9  2.731848e-11 4.431879e-14 4.431879e-14 4.431879e-14 4.431879e-14
#> 10 8.541631e-15 1.283609e-10 1.283609e-10 1.283609e-10 1.283609e-10
#> 11 9.272037e-14 5.402397e-12 5.402397e-12 5.402397e-12 5.402397e-12
#>                f
#> 1  -1.224145e-16
#> 2  -1.224145e-16
#> 3  -1.224145e-16
#> 4  -1.224145e-16
#> 5  -1.224145e-16
#> 6  -1.224145e-16
#> 7  -1.224145e-16
#> 8  -1.224145e-16
#> 9  -1.224145e-16
#> 10 -1.224145e-16
#> 11 -1.224145e-16

func1 <- function(x){ sin(10*x) - exp(-x) }
curve(func1,from=0,to=5)


x <- 2.04
numd1 <- nlmixr2Gill83(func1, x)
exact <- 10*cos(10*x) + exp(-x)
c(numd1$df, exact, (numd1$df - exact)/exact)
#> [1]  0.332398077  0.333537144 -0.003415112

x <- c(1:10)
numd1 <- nlmixr2Gill83(func1, x)
exact <- 10*cos(10*x) + exp(-x)
cbind(numd1=numd1$df, exact, err=(numd1$df - exact)/exact)
#>           numd1     exact           err
#>  [1,] -8.022836 -8.022836 -1.369260e-11
#>  [2,]  4.216156  4.216156 -2.839580e-11
#>  [3,]  1.592302  1.592302 -1.150871e-11
#>  [4,] -6.651065 -6.651065 -4.125002e-11
#>  [5,]  9.656398  9.656398 -5.172856e-11
#>  [6,] -9.521651 -9.521651  2.985064e-10
#>  [7,]  6.334104  6.334104 -8.697948e-11
#>  [8,] -1.103537 -1.103537 -9.731425e-11
#>  [9,] -4.480613 -4.480613 -1.320695e-10
#> [10,]  8.623852  8.623234  7.167430e-05

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
#> [1] 0.0010448
```
