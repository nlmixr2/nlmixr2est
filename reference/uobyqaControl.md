# Control for uobyqa estimation method in nlmixr2

Control for uobyqa estimation method in nlmixr2

## Usage

``` r
uobyqaControl(
  npt = NULL,
  rhobeg = NULL,
  rhoend = NULL,
  iprint = 0L,
  maxfun = 100000L,
  returnUobyqa = FALSE,
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
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  addProp = c("combined2", "combined1"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  eventSens = c("jump", "fd"),
  ...
)
```

## Arguments

- npt:

  Number of points for bobyqa's quadratic approximation to the
  objective; must be in \`\[n+2, (n+1)(n+2)/2\]\`. Defaults to \`2\*n +
  1\`. (bobyqa)

- rhobeg:

  Initial trust region radius for the bobyqa outer optimizer (with
  \`rhoend\`, must satisfy \`0 \< rhoend \< rhobeg\`). Default \`0.2\`
  (20 \`abs(upper-lower)/2\`. (bobyqa)

- rhoend:

  Final trust region radius. If not defined, \`10^(-sigdig-1)\` is used.
  (bobyqa)

- iprint:

  Controls amount of printing (\`0\`=none, \`1\`=start/end only,
  \`2\`=each new rho, \`3\`=every function evaluation, \`\>3\`=every
  \`iprint\` evaluations). Default \`0\`.

- maxfun:

  The maximum allowed number of function evaluations. If this is
  exceeded, the method will terminate.

- returnUobyqa:

  return the uobyqa output instead of the nlmixr2 fit

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

- normType:

  Parameter normalization/scaling used to get scaled initial values for
  `scaleType`, of the form `Vscaled = (Vunscaled-C1)/C2` (see [Feature
  Scaling](https://en.wikipedia.org/wiki/Feature_scaling); `rescale2`
  follows the
  [OptdesX](http://apmonitor.com/me575/uploads/Main/optimization_book.pdf)
  manual): `"rescale2"` scales all parameters to (-1, 1); `"rescale"`
  (min-max) scales to (0, 1); `"mean"` centers on the mean with range
  (0, 1); `"std"` standardizes by mean/sd; `"len"` scales to unit
  (Euclidean) length; `"constant"` performs no normalization (`C1=0`,
  `C2=1`).

- scaleType:

  The scaling scheme for nlmixr2: `"nlmixr2"` (default) scales as
  `(current-init)*scaleC[i] + scaleTo`, with `scaleTo` from `normType`
  and scales from `scaleC`; `"norm"` uses the simple scaling from
  `normType`; `"mult"` scales multiplicatively as
  `current/init*scaleTo`; `"multAdd"` scales linearly
  (`(current-init)+scaleTo`) for parameters in an exponential block
  (e.g. `exp(theta)`) and multiplicatively otherwise.

- scaleCmax:

  Maximum value of the scaleC to prevent overflow.

- scaleCmin:

  Minimum value of the scaleC to prevent underflow.

- scaleC:

  Scaling constant used with `scaleType="nlmixr2"`; when not specified,
  chosen by parameter type to keep gradient sizes similar on a log
  scale: \`1\` for exp()-transformed/power/boxCox/ yeoJohnson
  parameters, \`0.5\*abs(est)\` for additive/proportional/ lognormal
  error parameters, \`abs(1/digamma(est+1))\` for factorials, and
  \`log(abs(est))\*abs(est)\` for log-scale parameters. May be set
  explicitly per parameter if these defaults don't apply well.

- scaleTo:

  Scale the initial parameter estimate to this value. By default this
  is 1. When zero or below, no scaling is performed.

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

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- covMethod:

  Method for calculating the covariance. `"analytic"` (the default) uses
  the exact analytic observed-information R-matrix (reported as
  \\R^{-1}\\) and additionally returns the residual and `Omega` standard
  errors; it covers FOCEI/FOCE fits with additive, proportional, or
  combined error, mu-referenced/covariate/other structural parameters
  (and non-mu-referenced etas), and SD-scale inter-occasion variability,
  and emits a message and falls back to the finite-difference Hessian
  for anything out of scope (FO, `nAGQ > 1`, censoring, DV-transformed
  error, bounded-parameter transforms, a structural theta shared by two
  etas, non-SD `iovXform`, or a pure-proportional variance that vanishes
  at a near-zero prediction). The finite-difference methods use R (the
  Hessian) and S (the sum of individual gradient cross-products at the
  empirical Bayes estimates): `"r,s"` sandwich
  (`solve(R)%*%S%*%solve(R)`), `"r"` Hessian-based (`solve(R)`), `"s"`
  cross-product-based (`solve(S)`), or `""` to skip the covariance step.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- boundedTransform:

  When \`TRUE\` (default), bounded parameters are transformed for
  unbounded optimization methods and back-transformed for final
  estimates. \`FALSE\` optimizes on the original scale with bounds
  passed to the optimizer. \`NA\` transforms for optimization but skips
  the final back-transform.

- eventSens:

  Controls how dosing/event-parameter (\`alag\`, \`F\`, \`rate\`,
  \`dur\`) sensitivities are computed for THETA/ETA gradients:
  \`"jump"\` (default) uses rxode2's analytic event sensitivities;
  \`"fd"\` uses the legacy finite-difference behavior.

- ...:

  Ignored parameters

## Value

uobyqa control structure

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

fit2 <- nlmixr(mod, dsn, est="uobyqa")
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → pruning branches (`if`/`else`) of population log-likelihood model...
#> ✔ done
#> → loading llik model into symengine environment...
#> → finding duplicate expressions in population log-likelihood model...
#> → optimizing duplicate expressions in population log-likelihood model...
#> ✔ done
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
#> ── nlmixr² log-likelihood uobyqa ──
#> 
#>           OBJF     AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -679.7168 1164.16 1178.884      -579.0801        407.0929        67.95285
#> 
#> ── Time (sec $time): ──
#> 
#>             setup  optimize covariance preprocess postprocess table compress
#> elapsed 0.3307071 0.3253395   5.84e-06      0.045       0.013 0.024        0
#>             other
#> elapsed 0.1069476
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.    SE  %RSE  Back-transformed(95%CI) BSV(SD) Shrink(SD)%
#> E0  -0.8062 0.247 30.64 -0.8062 (-1.29, -0.3221)                    
#> Em    5.487 2.368 43.16    5.487 (0.8456, 10.13)                    
#> E50   2.898 1.244 42.94    2.898 (0.4592, 5.337)                    
#> g         2 FIXED FIXED                        2                    
#>  
#>   Covariance Type ($covMethod): r
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0383     0 -0.369 -0.805
#> 2 1     0.0509     0 -0.370 -0.805
#> 3 1     0.0543     1 -1.17  -0.804
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$nlm

fit2$uobyqa
#> parameter estimates: -0.806218115844621, 5.48662262692952, 2.89816168968597 
#> objective: 579.080149823443 
#> number of function evaluations: 66 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
