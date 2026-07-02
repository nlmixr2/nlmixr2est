# Control for lbfgsb3c estimation method in nlmixr2

Control for lbfgsb3c estimation method in nlmixr2

## Usage

``` r
lbfgsb3cControl(
  trace = 0,
  factr = 1e+07,
  pgtol = 0,
  abstol = 0,
  reltol = 0,
  lmm = 5L,
  maxit = 10000L,
  returnLbfgsb3c = FALSE,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  useColor = NULL,
  printNcol = NULL,
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
  addProp = c("combined2", "combined1"),
  eventSens = c("jump", "fd"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  ...
)
```

## Arguments

- trace:

  If positive, tracing information on the progress of the optimization
  is produced. Higher values may produce more tracing information: for
  method "L-BFGS-B" there are six levels of tracing. (To understand
  exactly what these do see the source code: higher levels give more
  detail.)

- factr:

  controls the convergence of the "L-BFGS-B" method. Convergence occurs
  when the reduction in the objective is within this factor of the
  machine tolerance. Default is 1e7, that is a tolerance of about 1e-8.

- pgtol:

  helps control the convergence of the "L-BFGS-B" method. It is a
  tolerance on the projected gradient in the current search direction.
  This defaults to zero, when the check is suppressed.

- abstol:

  helps control the convergence of the "L-BFGS-B" method. It is an
  absolute tolerance difference in x values. This defaults to zero, when
  the check is suppressed.

- reltol:

  helps control the convergence of the "L-BFGS-B" method. It is an
  relative tolerance difference in x values. This defaults to zero, when
  the check is suppressed.

- lmm:

  is an integer giving the number of BFGS updates retained in the
  "L-BFGS-B" method, It defaults to 5.

- maxit:

  maximum number of iterations.

- returnLbfgsb3c:

  return the lbfgsb3c output instead of the nlmixr2 fit

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

- useColor:

  Logical (or \`NULL\`) — whether to emit ANSI bold/color escapes in the
  iteration print. \`NULL\` (the default) defers to
  \[iterPrintControl()\]'s default (\[crayon::has_color()\]). Equivalent
  to \`print = iterPrintControl(useColor = ...)\`.

- printNcol:

  Integer (or \`NULL\`) — number of parameter columns emitted per row
  before wrapping to a continuation row. \`NULL\` (the default) defers
  to \[iterPrintControl()\]'s default (\`floor((getOption("width") - 23)
  / 12)\`, which fits an 80-column terminal). Equivalent to \`print =
  iterPrintControl(ncol = ...)\`.

- print:

  Either a scalar print-frequency (\`0\` = suppress iteration output;
  \`1\` (default) = print every parameter evaluation; \`N\` = print
  every Nth evaluation), OR a pre-built \[iterPrintControl()\] object
  bundling all iteration-print options (column wrap, header cadence,
  color, simple/three-row mode). The scalar form is equivalent to
  \`iterPrintControl(every = print, ncol = printNcol, useColor =
  useColor)\`.

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

- eventSens:

  How sensitivities of dosing/event parameters (absorption lag time,
  bioavailability, infusion rate and duration, etc.) are computed.
  \`"fd"\` uses the legacy finite differences. \`"jump"\` (the default)
  uses the analytic event ("jump") sensitivities provided by \`rxode2\`,
  which add accuracy and can speed up the gradient/Hessian by avoiding
  the extra finite-difference solves for these parameters.

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

  Ignored parameters

## Value

bobqya control structure

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

fit2 <- nlmixr(mod, dsn, est="lbfgsb3c")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
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

print(fit2)
#> ── nlmixr² log-likelihood lbfgsb3c ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -1529.377 314.5005 329.2238      -154.2503    1.921362e+19               1
#> 
#> ── Time (sec $time): ──
#> 
#>              setup    optimize preprocess postprocess table compress    other
#> elapsed 0.01620575 0.002775227      0.042       0.014 0.025    0.001 1.738019
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>       Est.        SE      %RSE        Back-transformed(95%CI) BSV(SD)
#> E0  0.1543 4.745e+04 3.076e+07 0.1543 (-9.301e+04, 9.301e+04)        
#> Em  -2.616 1.083e-05 0.0004139        -2.616 (-2.616, -2.616)        
#> E50  1.196 1.348e-05  0.001128           1.196 (1.196, 1.196)        
#> g        2     FIXED     FIXED                              2        
#>     Shrink(SD)%
#> E0             
#> Em             
#> E50            
#> g              
#>  
#>   Covariance Type ($covMethod): r
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     ABNORMAL_TERMINATION_IN_LNSRCH 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED     v
#>   <fct>  <dbl> <dbl>  <dbl> <dbl>
#> 1 1     0.0399     1 -0.620 0.151
#> 2 1     0.0448     0 -0.771 0.151
#> 3 1     0.0681     0 -0.769 0.146
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$lbfgsb3c

fit2$lbfgsb3c
#> $par
#>         E0         Em        E50 
#>  0.1542503 -2.6155545  1.1955657 
#> 
#> $grad
#> [1] -1.3077772  0.3985219 -0.4258109
#> 
#> $value
#> [1] 154.2503
#> 
#> $counts
#> [1] 48 48
#> 
#> $convergence
#> [1] NA
#> 
#> $message
#> [1] "ABNORMAL_TERMINATION_IN_LNSRCH"
#> 
#> $scaleC
#> [1] 0.0005000000 0.0003333333 0.0004151442
#> 
#> $par.scaled
#>         E0         Em        E50 
#>  -692.4995 -9347.6634 -1936.7227 
#> 
#> $hessian
#>           E0 Em E50
#> E0  -3792.25  0   0
#> Em      0.00  0   0
#> E50     0.00  0   0
#> 
#> $cov.scaled
#>               E0          Em         E50
#> E0  9.007199e+15 0.000000000 0.000000000
#> Em  0.000000e+00 0.001054783 0.000000000
#> E50 0.000000e+00 0.000000000 0.001054783
#> 
#> $r
#>            E0 Em E50
#> E0  -1896.125  0   0
#> Em      0.000  0   0
#> E50     0.000  0   0
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
