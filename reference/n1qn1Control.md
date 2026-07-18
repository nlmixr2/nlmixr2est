# Control for n1qn1 estimation method in nlmixr2

Control for n1qn1 estimation method in nlmixr2

## Usage

``` r
n1qn1Control(
  epsilon = (.Machine$double.eps)^0.25,
  max_iterations = 10000,
  nsim = 10000,
  imp = 0,
  print.functions = FALSE,
  returnN1qn1 = FALSE,
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
  sensMethod = c("default", "forward", "adjoint"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", "n1qn1", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  ...
)
```

## Arguments

- epsilon:

  Precision of estimate for n1qn1 optimization.

- max_iterations:

  Number of iterations

- nsim:

  Number of function evaluations

- imp:

  Verbosity of messages.

- print.functions:

  Boolean to control if the function value and parameter estimates are
  echoed every time a function is called.

- returnN1qn1:

  return the n1qn1 output instead of the nlmixr2 fit

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

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- eventSens:

  Controls how dosing/event-parameter (\`alag\`, \`F\`, \`rate\`,
  \`dur\`) sensitivities are computed for THETA/ETA gradients:
  \`"jump"\` (default) uses rxode2's analytic event sensitivities;
  \`"fd"\` uses the legacy finite-difference behavior.

- sensMethod:

  Method used to compute the ODE parameter sensitivities: \`"default"\`
  (the default) defers to the global option
  \`getOption("nlmixr2est.adjoint")\` (itself \`"forward"\` by default);
  \`"forward"\` uses the classic variational (forward) sensitivity ODEs;
  \`"adjoint"\` uses the in-engine discrete adjoint with the matching
  adjoint (\`s\`) method.

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

fit2 <- nlmixr(mod, dsn, est="n1qn1")
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
#> ── nlmixr² log-likelihood n1qn1 ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -715.2361 1128.641 1143.364      -561.3205        493.1807        71.78878
#> 
#> ── Time (sec $time): ──
#> 
#>             setup  optimize covariance preprocess postprocess table compress
#> elapsed 0.3388906 0.5453451  6.722e-06      0.053       0.006 0.026    0.001
#>             other
#> elapsed 0.1057575
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI)
#> E0  -0.7540 0.2468 32.73 -0.7540 (-1.238, -0.2703)
#> Em    6.000  2.701 45.01     6.000 (0.7071, 11.29)
#> E50   2.931  1.257 42.87     2.931 (0.4682, 5.395)
#> g     2.000  FIXED FIXED                     2.000
#>  
#>   Covariance Type ($covMethod): r
#>   Censoring ($censInformation): No censoring
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0591     0 -0.386 -0.752
#> 2 1     0.0591     0 -0.386 -0.752
#> 3 1     0.0628     0 -0.386 -0.751
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$n1qn1

fit2$n1qn1
#> $value
#> [1] 561.3205
#> 
#> $par
#>        E0        Em       E50 
#> -0.753981  6.000383  2.931381 
#> 
#> $H
#>              [,1]         [,2]         [,3]
#> [1,]  0.001650649  0.002660560 -0.007137467
#> [2,]  0.002660560  0.009466189 -0.020792349
#> [3,] -0.007137467 -0.020792349  0.050431465
#> 
#> $c.hess
#>  [1]  0.001650649  0.002660560 -0.007137467  0.009466189 -0.020792349
#>  [6]  0.050431465  0.000000000  0.000000000  0.000000000  0.000000000
#> [11]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#> [16]  0.000000000  0.000000000  0.000000000  0.000000000  0.000000000
#> [21]  0.000000000  0.000000000  0.000000000  0.000000000
#> 
#> $n.fn
#> [1] 37
#> 
#> $n.gr
#> [1] 37
#> 
#> $scaleC
#> [1] 0.002875081 0.036929293 0.033369905
#> 
#> $par.scaled
#>        E0        Em       E50 
#> -437.1549  147.9436   28.9108 
#> 
#> $hessian
#>               E0           Em          E50
#> E0   0.001602122  0.002609465 -0.006942853
#> Em   0.002609465  0.009410840 -0.020574667
#> E50 -0.006942853 -0.020574667  0.049546009
#> 
#> $cov.scaled
#>           E0       Em      E50
#> E0  7368.597 2325.947 1998.439
#> Em  2325.947 5348.213 2546.853
#> E50 1998.439 2546.853 1418.389
#> 
#> $r
#>                E0           Em          E50
#> E0   0.0008010608  0.001304733 -0.003471426
#> Em   0.0013047326  0.004705420 -0.010287333
#> E50 -0.0034714263 -0.010287333  0.024773004
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
