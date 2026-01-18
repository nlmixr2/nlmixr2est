# nlmixr2 optim defaults

nlmixr2 optim defaults

## Usage

``` r
optimControl(
  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
  trace = 0,
  fnscale = 1,
  parscale = 1,
  ndeps = 0.001,
  maxit = 10000,
  abstol = 1e-08,
  reltol = 1e-08,
  alpha = 1,
  beta = 0.5,
  gamma = 2,
  REPORT = NULL,
  warn.1d.NelderMead = TRUE,
  type = NULL,
  lmm = 5,
  factr = 1e+07,
  pgtol = 0,
  temp = 10,
  tmax = 10,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1/3),
  shi21maxFD = 20L,
  solveType = c("grad", "fun"),
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
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  returnOptim = FALSE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", "optim", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  ...
)
```

## Arguments

- method:

  The method to be used. See ‘Details’. Can be abbreviated.

- trace:

  Non-negative integer. If positive, tracing information on the progress
  of the optimization is produced. Higher values may produce more
  tracing information: for method \`"L-BFGS-B"\`, there are six levels
  of tracing. See \`optim()\` for more information

- fnscale:

  An overall scaling to be applied to the value of \`fn\` and \`gr\`
  during optimization. If negative, turns the problem into a
  maximization problem. Optimization is performed on \`fn(par)/fnscale\`

- parscale:

  A vector of scaling values for the parameters. Optimization is
  performed on \`par/parscale\` and these should be comparable in the
  sense that a unit change in any element produces about a unit change
  in the scaled value. Not used (nor needed) for \`method = "Brent"\`

- ndeps:

  A vector of step sizes for the finite-difference approximation to the
  gradient, on \`par/parscale\` scale. Defaults to \`1e-3\`

- maxit:

  The maximum number of iterations. Defaults to \`100\` for the
  derivative-based methods, and \`500\` for \`"Nelder-Mead"\`.

- abstol:

  The absolute convergence tolerance. Only useful for non-negative
  functions, as a tolerance for reaching zero.

- reltol:

  Relative convergence tolerance. The algorithm stops if it is unable to
  reduce the value by a factor of \`reltol \* (abs(val) + reltol)\` at a
  step

- alpha:

  Reflection factor for the \`"Nelder-Mead"\` method.

- beta:

  Contraction factor for the \`"Nelder-Mead"\` method

- gamma:

  Expansion factor for the \`"Nelder-Mead"\` method

- REPORT:

  The frequency of reports for the \`"BFGS"\`, \`"L-BFGS-B"\` and
  \`"SANN"\` methods if \`control\$trace\` is positive. Defaults to
  every 10 iterations for \`"BFGS"\` and \`"L-BFGS-B"\`, or every 100
  temperatures for \`"SANN"\`

- warn.1d.NelderMead:

  a logical indicating if the (default) \`"Nelder-Mead"\` method should
  signal a warning when used for one-dimensional minimization. As the
  warning is sometimes inappropriate, you can suppress it by setting
  this option to \`FALSE\`

- type:

  for the conjugate-gradients method. Takes value \`1\` for the
  Fletcher-Reeves update, \`2\` for Polak-Ribiere and \`3\` for
  Beale-Sorenson.

- lmm:

  is an integer giving the number of BFGS updates retained in the
  \`"L-BFGS-B"\` method, It defaults to \`5\`

- factr:

  controls the convergence of the \`"L-BFGS-B"\` method. Convergence
  occurs when the reduction in the objective is within this factor of
  the machine tolerance. Default is \`1e7\`, that is a tolerance of
  about \`1e-8\`.

- pgtol:

  helps control the convergence of the ‘"L-BFGS-B"’ method. It is a
  tolerance on the projected gradient in the current search direction.
  This defaults to zero, when the check is suppressed

- temp:

  controls the \`"SANN"\` method. It is the starting temperature for the
  cooling schedule. Defaults to \`10\`.

- tmax:

  is the number of function evaluations at each temperature for the
  \`"SANN"\` method. Defaults to \`10\`.

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

- solveType:

  tells if \`optim\` will use nlmixr2's analytical gradients when
  available (finite differences will be used for event-related
  parameters like parameters controlling lag time, duration/rate of
  infusion, and modeled bioavailability). This can be:

  \- \`"gradient"\` which will use the gradient and let \`optim\`
  calculate the finite difference hessian

  \- \`"fun"\` where optim will calculate both the finite difference
  gradient and the finite difference Hessian

  When using nlmixr2's finite differences, the "ideal" step size for
  either central or forward differences are optimized for with the
  Shi2021 method which may give more accurate derivatives

  These are only applied in the gradient based methods: "BFGS", "CG",
  "L-BFGS-B"

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

- returnOptim:

  logical; when TRUE this will return the optim list instead of the
  nlmixr2 fit object

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

  allows selection of "r", which uses nlmixr2's \`nlmixr2Hess()\` for
  the hessian calculation or "optim" which uses the hessian from
  \`stats::optim(.., hessian=TRUE)\`

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

  Further arguments to be passed to `fn` and `gr`.

## Value

optimControl object for nlmixr2

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

fit2 <- nlmixr(mod, dsn, est="optim", optimControl(method="BFGS"))
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
fit2
#> ── nlmixr² log-likelihood optim with BFGS method ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -701.3632 1142.514 1157.237      -568.2569        253.4253        43.23195
#> 
#> ── Time (sec value$time): ──
#> 
#>            setup table compress    other
#> elapsed 0.003121 0.044    0.001 2.765879
#> 
#> ── (value$parFixed or value$parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> E0  -0.7564 0.2336 30.88 -0.7564 (-1.214, -0.2987)                    
#> Em    5.329  1.895 35.56      5.329 (1.615, 9.044)                    
#> E50   2.745 0.9886 36.01     2.745 (0.8076, 4.683)                    
#> g         2  FIXED FIXED                         2                    
#>  
#>   Covariance Type (value$covMethod): r (optim)
#>   Censoring (value$censInformation): No censoring
#> 
#> ── Fit Data (object value is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0429     0 -0.385 -0.755
#> 2 1     0.0693     0 -0.386 -0.753
#> 3 1     0.0743     0 -0.386 -0.753
#> # ℹ 997 more rows
# }
```
