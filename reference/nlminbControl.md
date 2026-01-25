# nlmixr2 nlminb defaults

nlmixr2 nlminb defaults

## Usage

``` r
nlminbControl(
  eval.max = 200,
  iter.max = 150,
  trace = 0,
  abs.tol = 0,
  rel.tol = 1e-10,
  x.tol = 1.5e-08,
  xf.tol = 2.2e-14,
  step.min = 1,
  step.max = 1,
  sing.tol = rel.tol,
  scale = 1,
  scale.init = NULL,
  diff.g = NULL,
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  returnNlminb = FALSE,
  solveType = c("hessian", "grad", "fun"),
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1/3),
  shi21maxFD = 20L,
  optimHessType = c("central", "forward"),
  hessErr = (.Machine$double.eps)^(1/3),
  shi21maxHess = 20L,
  useColor = crayon::has_color(),
  printNcol = floor((getOption("width") - 23)/12),
  print = 1L,
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
  scaleType = c("nlmixr2", "norm", "mult", "multAdd"),
  scaleCmax = 1e+05,
  scaleCmin = 1e-05,
  scaleC = NULL,
  scaleTo = 1,
  gradTo = 1,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = TRUE,
  covMethod = c("r", "nlminb", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  ...
)
```

## Arguments

- eval.max:

  Maximum number of evaluations of the objective function allowed.
  Defaults to 200.

- iter.max:

  Maximum number of iterations allowed. Defaults to 150.

- trace:

  The value of the objective function and the parameters is printed
  every trace'th iteration. When 0 no trace information is to be printed

- abs.tol:

  Absolute tolerance. Defaults to 0 so the absolute convergence test is
  not used. If the objective function is known to be non-negative, the
  previous default of \`1e-20\` would be more appropriate

- rel.tol:

  Relative tolerance. Defaults to \`1e-10\`.

- x.tol:

  X tolerance. Defaults to \`1.5e-8\`.

- xf.tol:

  false convergence tolerance. Defaults to \`2.2e-14\`.

- step.min:

  Minimum step size. Default to ‘1.’.

- step.max:

  Maximum step size. Default to ‘1.’.

- sing.tol:

  singular convergence tolerance; defaults to \`rel.tol;.

- scale:

  See PORT documentation (or leave alone).

- scale.init:

  ... probably need to check PORT documentation

- diff.g:

  an estimated bound on the relative error in the objective function
  value

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- returnNlminb:

  logical; when TRUE this will return the nlminb result instead of the
  nlmixr2 fit object

- solveType:

  tells if \`nlm\` will use nlmixr2's analytical gradients when
  available (finite differences will be used for event-related
  parameters like parameters controlling lag time, duration/rate of
  infusion, and modeled bioavailability). This can be:

  \- \`"hessian"\` which will use the analytical gradients to create a
  Hessian with finite differences.

  \- \`"gradient"\` which will use the gradient and let \`nlm\`
  calculate the finite difference hessian

  \- \`"fun"\` where nlm will calculate both the finite difference
  gradient and the finite difference Hessian

  When using nlmixr2's finite differences, the "ideal" step size for
  either central or forward differences are optimized for with the
  Shi2021 method which may give more accurate derivatives

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- eventType:

  Event gradient type for dosing events; Can be "central" or "forward"

- shiErr:

  This represents the epsilon when optimizing the ideal step size for
  numeric differentiation using the Shi2021 method

- shi21maxFD:

  The maximum number of steps for the optimization of the forward
  difference step size when using dosing events (lag time, modeled
  duration/rate and bioavailability)

- optimHessType:

  The hessian type for when calculating the individual hessian by
  numeric differences (in generalized log-likelihood estimation). The
  options are "central", and "forward". The central differences is what
  R's \`optimHess()\` uses and is the default for this method. (Though
  the "forward" is faster and still reasonable for most cases). The
  Shi21 cannot be changed for the Gill83 algorithm with the optimHess in
  a generalized likelihood problem.

- hessErr:

  This represents the epsilon when optimizing the Hessian step size
  using the Shi2021 method.

- shi21maxHess:

  Maximum number of times to optimize the best step size for the hessian
  calculation

- useColor:

  Boolean indicating if focei can use ASCII color codes

- printNcol:

  Number of columns to printout before wrapping parameter
  estimates/gradient

- print:

  Integer representing when the outer step is printed. When this is 0 or
  do not print the iterations. 1 is print every function evaluation
  (default), 5 is print every 5 evaluations.

- normType:

  This is the type of parameter normalization/scaling used to get the
  scaled initial values for nlmixr2. These are used with `scaleType` of.

  With the exception of `rescale2`, these come from [Feature
  Scaling](https://en.wikipedia.org/wiki/Feature_scaling). The
  `rescale2` The rescaling is the same type described in the
  [OptdesX](http://apmonitor.com/me575/uploads/Main/optimization_book.pdf)
  software manual.

  In general, all all scaling formula can be described by:

  \$\$v\_{scaled}\$\$ = (\$\$v\_{unscaled}-C\_{1}\$\$)/\$\$C\_{2}\$\$

  Where

  The other data normalization approaches follow the following formula

  \$\$v\_{scaled}\$\$ = (\$\$v\_{unscaled}-C\_{1}\$\$)/\$\$C\_{2}\$\$

  - `rescale2` This scales all parameters from (-1 to 1). The relative
    differences between the parameters are preserved with this approach
    and the constants are:

    \$\$C\_{1}\$\$ = (max(all unscaled values)+min(all unscaled
    values))/2

    \$\$C\_{2}\$\$ = (max(all unscaled values) - min(all unscaled
    values))/2

  - `rescale` or min-max normalization. This rescales all parameters
    from (0 to 1). As in the `rescale2` the relative differences are
    preserved. In this approach:

    \$\$C\_{1}\$\$ = min(all unscaled values)

    \$\$C\_{2}\$\$ = max(all unscaled values) - min(all unscaled values)

  - `mean` or mean normalization. This rescales to center the parameters
    around the mean but the parameters are from 0 to 1. In this
    approach:

    \$\$C\_{1}\$\$ = mean(all unscaled values)

    \$\$C\_{2}\$\$ = max(all unscaled values) - min(all unscaled values)

  - `std` or standardization. This standardizes by the mean and standard
    deviation. In this approach:

    \$\$C\_{1}\$\$ = mean(all unscaled values)

    \$\$C\_{2}\$\$ = sd(all unscaled values)

  - `len` or unit length scaling. This scales the parameters to the unit
    length. For this approach we use the Euclidean length, that is:

    \$\$C\_{1}\$\$ = 0

    \$\$C\_{2}\$\$ = \$\$\sqrt(v_1^2 + v_2^2 + \cdots + v_n^2)\$\$

  - `constant` which does not perform data normalization. That is

    \$\$C\_{1}\$\$ = 0

    \$\$C\_{2}\$\$ = 1

- scaleType:

  The scaling scheme for nlmixr2. The supported types are:

  - `nlmixr2` In this approach the scaling is performed by the following
    equation:

    \$\$v\_{scaled}\$\$ = (\$\$v\_{current} -
    v\_{init}\$\$)\*scaleC\[i\] + scaleTo

    The `scaleTo` parameter is specified by the `normType`, and the
    scales are specified by `scaleC`.

  - `norm` This approach uses the simple scaling provided by the
    `normType` argument.

  - `mult` This approach does not use the data normalization provided by
    `normType`, but rather uses multiplicative scaling to a constant
    provided by the `scaleTo` argument.

    In this case:

    \$\$v\_{scaled}\$\$ =
    \$\$v\_{current}\$\$/\$\$v\_{init}\$\$\*scaleTo

  - `multAdd` This approach changes the scaling based on the parameter
    being specified. If a parameter is defined in an exponential block
    (ie exp(theta)), then it is scaled on a linearly, that is:

    \$\$v\_{scaled}\$\$ = (\$\$v\_{current}-v\_{init}\$\$) + scaleTo

    Otherwise the parameter is scaled multiplicatively.

    \$\$v\_{scaled}\$\$ =
    \$\$v\_{current}\$\$/\$\$v\_{init}\$\$\*scaleTo

- scaleCmax:

  Maximum value of the scaleC to prevent overflow.

- scaleCmin:

  Minimum value of the scaleC to prevent underflow.

- scaleC:

  The scaling constant used with `scaleType=nlmixr2`. When not
  specified, it is based on the type of parameter that is estimated. The
  idea is to keep the derivatives similar on a log scale to have similar
  gradient sizes. Hence parameters like log(exp(theta)) would have a
  scaling factor of 1 and log(theta) would have a scaling factor of
  ini_value (to scale by 1/value; ie d/dt(log(ini_value)) = 1/ini_value
  or scaleC=ini_value)

  - For parameters in an exponential (ie exp(theta)) or parameters
    specifying powers, boxCox or yeoJohnson transformations , this is 1.

  - For additive, proportional, lognormal error structures, these are
    given by 0.5\*abs(initial_estimate)

  - Factorials are scaled by abs(1/digamma(initial_estimate+1))

  - parameters in a log scale (ie log(theta)) are transformed by
    log(abs(initial_estimate))\*abs(initial_estimate)

  These parameter scaling coefficients are chose to try to keep similar
  slopes among parameters. That is they all follow the slopes
  approximately on a log-scale.

  While these are chosen in a logical manner, they may not always apply.
  You can specify each parameters scaling factor by this parameter if
  you wish.

- scaleTo:

  Scale the initial parameter estimate to this value. By default this
  is 1. When zero or below, no scaling is performed.

- gradTo:

  this is the factor that the gradient is scaled to before optimizing.
  This only works with scaleType="nlmixr2".

- addProp:

  specifies the type of additive plus proportional errors, the one where
  standard deviations add (combined1) or the type where the variances
  add (combined2).

  The combined1 error type can be described by the following equation:

  \$\$y = f + (a + b\times f^c) \times \varepsilon\$\$

  The combined2 error model can be described by the following equation:

  \$\$y = f + \sqrt{a^2 + b^2\times f^{2\times c}} \times
  \varepsilon\$\$

  Where:

  \- y represents the observed value

  \- f represents the predicted value

  \- a is the additive standard deviation

  \- b is the proportional/power standard deviation

  \- c is the power exponent (in the proportional case c=1)

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- covMethod:

  Method for calculating covariance. In this discussion, R is the
  Hessian matrix of the objective function. The S matrix is the sum of
  individual gradient cross-product (evaluated at the individual
  empirical Bayes estimates).

  - "`r,s`" Uses the sandwich matrix to calculate the covariance, that
    is: `solve(R) %*% S %*% solve(R)`

  - "`r`" Uses the Hessian matrix to calculate the covariance as
    `2 %*% solve(R)`

  - "`s`" Uses the cross-product matrix to calculate the covariance as
    `4 %*% solve(S)`

  - "" Does not calculate the covariance step.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits. This controls:

  - The tolerance of the inner and outer optimization is `10^-sigdig`

  - The tolerance of the ODE solvers is `0.5*10^(-sigdig-2)`; For the
    sensitivity equations and steady-state solutions the default is
    `0.5*10^(-sigdig-1.5)` (sensitivity changes only applicable for
    liblsoda)

  - The tolerance of the boundary check is `5 * 10 ^ (-sigdig + 1)`

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- ...:

  Further arguments to be supplied to `objective`.

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{
# A logit regression example with emax model

dsn <- data.frame(i=1:1000)
dsn$time <- exp(rnorm(1000))
dsn$DV=rbinom(1000,1,exp(-1+dsn$time)/(1+exp(-1+dsn$time)))

mod <- function() {
 ini({
   E0 <- 0.5
   Em <- 0.5
   E50 <- 2
   g <- fix(2)
 })
 model({
   v <- E0+Em*time^g/(E50^g+time^g)
   ll(bin) ~ DV * v - log(1 + exp(v))
 })
}

fit2 <- nlmixr(mod, dsn, est="nlminb")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
#> → calculate jacobian
#> → calculate ∂(f)/∂(θ)
#> → finding duplicate expressions in nlm llik gradient...
#> → optimizing duplicate expressions in nlm llik gradient...
#> → finding duplicate expressions in nlm pred-only...
#> → optimizing duplicate expressions in nlm pred-only...
#>  
#>  
#>  
#>  
#> → calculating covariance
#> ✔ done
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of full model...
#> ✔ done
#> → finding duplicate expressions in EBE model...
#> → optimizing duplicate expressions in EBE model...
#> → compiling EBE model...
#>  
#>  
#> ✔ done
#> → Calculating residuals/tables
#> ✔ done
#> → compress origData in nlmixr2 object, save 8304
#> → compress parHistData in nlmixr2 object, save 2704

print(fit2)
#> ── nlmixr² log-likelihood nlminb ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -666.1358 1177.741 1192.465      -585.8706        1303.943        97.41604
#> 
#> ── Time (sec $time): ──
#> 
#>            setup table compress    other
#> elapsed 0.002983 0.039    0.012 1.100017
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE    Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> E0  -0.5491 0.2149 39.14 -0.5491 (-0.9703, -0.1279)                    
#> Em    6.749  4.472 66.27      6.749 (-2.017, 15.51)                    
#> E50   3.802  2.055 54.04     3.802 (-0.2252, 7.829)                    
#> g         2  FIXED FIXED                          2                    
#>  
#>   Covariance Type ($covMethod): r (nlminb)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     relative convergence (4) 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0568     1 -1.00  -0.548
#> 2 1     0.0638     0 -0.457 -0.547
#> 3 1     0.0685     1 -1.00  -0.547
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$nlminb

fit2$nlminb
#> $par
#>         E0         Em        E50 
#> -0.5491283  6.7490353  3.8019249 
#> 
#> $objective
#> [1] 585.8706
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 9
#> 
#> $evaluations
#> function gradient 
#>       16       10 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
#> $scaleC
#> [1] 0.002838234 0.034376196 0.031852867
#> 
#> $parHistData
#>    iter                type     objf            E0            Em           E50
#> 1     1              Scaled 680.3819 -1.000000e+00 -1.000000e+00  1.000000e+00
#> 2     1            Unscaled 680.3819  5.000000e-01  5.000000e-01  2.000000e+00
#> 3     1    Back-Transformed 680.3819  5.000000e-01  5.000000e-01  2.000000e+00
#> 4     2              Scaled 679.6912 -1.356390e+00 -1.487991e-01  1.098202e+00
#> 5     2            Unscaled 679.6912  4.989885e-01  5.292610e-01  2.003128e+00
#> 6     2    Back-Transformed 679.6912  4.989885e-01  5.292610e-01  2.003128e+00
#> 7     3              Scaled 677.7292 -2.662247e+00  2.312230e+00  1.559640e+00
#> 8     3            Unscaled 677.7292  4.952822e-01  6.138619e-01  2.017826e+00
#> 9     3    Back-Transformed 677.7292  4.952822e-01  6.138619e-01  2.017826e+00
#> 10    4              Scaled 671.1568 -1.052566e+01  9.853269e+00  5.228357e+00
#> 11    4            Unscaled 671.1568  4.729640e-01  8.730941e-01  2.134685e+00
#> 12    4    Back-Transformed 671.1568  4.729640e-01  8.730941e-01  2.134685e+00
#> 13    5              Scaled 651.0013 -5.190061e+01  2.540722e+01  2.199849e+01
#> 14    5            Unscaled 651.0013  3.555322e-01  1.407780e+00  2.668862e+00
#> 15    5    Back-Transformed 651.0013  3.555322e-01  1.407780e+00  2.668862e+00
#> 16    6              Scaled 614.1101 -2.075764e+02  6.034518e+01  5.593261e+01
#> 17    6            Unscaled 614.1101 -8.631220e-02  2.608814e+00  3.749761e+00
#> 18    6    Back-Transformed 614.1101 -8.631220e-02  2.608814e+00  3.749761e+00
#> 19    7              Scaled 738.6608 -3.752646e+02  9.912716e+01 -9.016151e+01
#> 20    7            Unscaled 738.6608 -5.622504e-01  3.941991e+00 -9.037553e-01
#> 21    7    Back-Transformed 738.6608 -5.622504e-01  3.941991e+00 -9.037553e-01
#> 22    8              Scaled 600.5570 -2.261833e+02  8.595275e+01  2.977248e+01
#> 23    8            Unscaled 600.5570 -1.391227e-01  3.489105e+00  2.916486e+00
#> 24    8    Back-Transformed 600.5570 -1.391227e-01  3.489105e+00  2.916486e+00
#> 25    9              Scaled 593.9449 -2.615096e+02  1.052749e+02  3.988392e+01
#> 26    9            Unscaled 593.9449 -2.393872e-01  4.153327e+00  3.238564e+00
#> 27    9    Back-Transformed 593.9449 -2.393872e-01  4.153327e+00  3.238564e+00
#> 28   10              Scaled 586.7110 -3.573225e+02  1.362035e+02  4.340594e+01
#> 29   10            Unscaled 586.7110 -5.113266e-01  5.216534e+00  3.350751e+00
#> 30   10    Back-Transformed 586.7110 -5.113266e-01  5.216534e+00  3.350751e+00
#> 31   11              Scaled 585.9822 -3.789858e+02  1.547439e+02  4.437743e+01
#> 32   11            Unscaled 585.9822 -5.728119e-01  5.853884e+00  3.381696e+00
#> 33   11    Back-Transformed 585.9822 -5.728119e-01  5.853884e+00  3.381696e+00
#> 34   12              Scaled 585.8810 -3.724411e+02  1.719870e+02  5.352939e+01
#> 35   12            Unscaled 585.8810 -5.542366e-01  6.446635e+00  3.673212e+00
#> 36   12    Back-Transformed 585.8810 -5.542366e-01  6.446635e+00  3.673212e+00
#> 37   13              Scaled 585.8709 -3.709141e+02  1.794151e+02  5.695939e+01
#> 38   13            Unscaled 585.8709 -5.499028e-01  6.701984e+00  3.782467e+00
#> 39   13    Back-Transformed 585.8709 -5.499028e-01  6.701984e+00  3.782467e+00
#> 40   14              Scaled 585.8706 -3.706539e+02  1.807306e+02  5.754553e+01
#> 41   14            Unscaled 585.8706 -5.491642e-01  6.747207e+00  3.801137e+00
#> 42   14    Back-Transformed 585.8706 -5.491642e-01  6.747207e+00  3.801137e+00
#> 43   15              Scaled 585.8706 -3.706415e+02  1.807830e+02  5.756985e+01
#> 44   15            Unscaled 585.8706 -5.491290e-01  6.749008e+00  3.801912e+00
#> 45   15    Back-Transformed 585.8706 -5.491290e-01  6.749008e+00  3.801912e+00
#> 46   16              Scaled 585.8706 -3.706413e+02  1.807838e+02  5.757026e+01
#> 47   16            Unscaled 585.8706 -5.491283e-01  6.749035e+00  3.801925e+00
#> 48   16    Back-Transformed 585.8706 -5.491283e-01  6.749035e+00  3.801925e+00
#> 49   17              Scaled 585.8706 -3.706413e+02  1.807838e+02  5.757026e+01
#> 50   17            Unscaled 585.8706 -5.491283e-01  6.749035e+00  3.801925e+00
#> 51   17    Back-Transformed 585.8706 -5.491283e-01  6.749035e+00  3.801925e+00
#> 52    1 Forward Sensitivity       NA  2.763308e-01 -7.019160e-01 -6.944356e-02
#> 53    6 Forward Sensitivity       NA  3.428367e-02 -5.252987e-01  3.019175e-01
#> 54    8 Forward Sensitivity       NA  1.571286e-01 -4.764074e-02 -1.376144e-01
#> 55   10 Forward Sensitivity       NA -2.448076e-04 -7.897922e-02  9.805760e-02
#> 56   11 Forward Sensitivity       NA -1.440626e-03 -9.054192e-03 -8.703912e-04
#> 57   12 Forward Sensitivity       NA -4.392492e-04 -3.885017e-03  3.221958e-03
#> 58   13 Forward Sensitivity       NA -1.138820e-04 -6.531890e-04  7.308369e-04
#> 59   14 Forward Sensitivity       NA -3.409368e-06 -1.849928e-05  1.409273e-05
#> 60   15 Forward Sensitivity       NA  2.857787e-09  1.857756e-08 -4.369005e-07
#> 61   16 Forward Sensitivity       NA  1.351871e-10  8.014625e-10 -8.162520e-09
#> 
#> $par.scaled
#>         E0         Em        E50 
#> -370.64126  180.78379   57.57026 
#> 
#> $hessian
#>               E0           Em          E50
#> E0   0.001648690  0.001674884 -0.004309052
#> Em   0.001674884  0.004421178 -0.009619118
#> E50 -0.004309052 -0.009619118  0.022325585
#> 
#> $covMethod
#> [1] "r (nlminb)"
#> 
#> $cov.scaled
#>           E0        Em      E50
#> E0  5733.154  3764.455 2728.493
#> Em  3764.455 16926.772 8019.581
#> E50 2728.493  8019.581 4161.079
#> 
#> $cov
#>             E0         Em       E50
#> E0  0.04618383  0.3672891 0.2466718
#> Em  0.36728908 20.0027536 8.7812844
#> E50 0.24667177  8.7812844 4.2218517
#> 
#> $r
#>                E0            Em          E50
#> E0   0.0008243449  0.0008374419 -0.002154526
#> Em   0.0008374419  0.0022105890 -0.004809559
#> E50 -0.0021545261 -0.0048095589  0.011162792
#> 
# }
```
