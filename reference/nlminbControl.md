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
  indTolRelax = TRUE,
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

- indTolRelax:

  When \`TRUE\` (default), only subjects whose ODE solve produced
  NaN/Inf have their tolerances relaxed, and the relaxed tolerance
  persists across optimizer calls (sticky). When \`FALSE\`, all subjects
  have their tolerances relaxed on each retry and tolerances are reset
  afterward.

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
#> → compress origData in nlmixr2 object, save 8360
#> → compress parHistData in nlmixr2 object, save 2904

print(fit2)
#> ── nlmixr² log-likelihood nlminb ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -634.9123 1208.965 1223.688      -601.4824          120009        935.5517
#> 
#> ── Time (sec $time): ──
#> 
#>              setup    optimize preprocess postprocess table compress     other
#> elapsed 0.01578497 0.000529957      0.046       0.012 0.023    0.012 0.8526851
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE     Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> E0  -0.3711 0.1918  51.7 -0.3711 (-0.7471, 0.004905)                    
#> Em    17.26  45.73   265       17.26 (-72.38, 106.9)                    
#> E50   7.909  12.21 154.3       7.909 (-16.01, 31.83)                    
#> g         2  FIXED FIXED                           2                    
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
#> 1 1     0.0261     0 -0.525 -0.371
#> 2 1     0.0476     0 -0.525 -0.370
#> 3 1     0.0505     0 -0.525 -0.370
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$nlminb

fit2$nlminb
#> $par
#>         E0         Em        E50 
#> -0.3710903 17.2550218  7.9093700 
#> 
#> $objective
#> [1] 601.4824
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 14
#> 
#> $evaluations
#> function gradient 
#>       21       15 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
#> $scaleC
#> [1] 0.002863784 0.030813757 0.029116520
#> 
#> $par.scaled
#>        E0        Em       E50 
#> -305.1746  542.7513  203.9559 
#> 
#> $hessian
#>                E0            Em          E50
#> E0   0.0017400253  0.0004542875 -0.001736241
#> Em   0.0004542875  0.0003568154 -0.001293392
#> E50 -0.0017362413 -0.0012933920  0.004717977
#> 
#> $cov.scaled
#>            E0         Em       E50
#> E0   4487.334   43384.11  13544.74
#> Em  43384.107 2202580.85 619783.89
#> E50 13544.738  619783.89 175740.70
#> 
#> $r
#>                E0            Em           E50
#> E0   0.0008700127  0.0002271437 -0.0008681206
#> Em   0.0002271437  0.0001784077 -0.0006466960
#> E50 -0.0008681206 -0.0006466960  0.0023589883
#> 
# }
```
