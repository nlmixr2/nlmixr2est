# Control for lbfgsb3c estimation method in nlmixr2

Control for lbfgsb3c estimation method in nlmixr2

## Usage

``` r
lbfgsb3cControl(
  trace = 0,
  factr = NULL,
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
  sensMethod = c("default", "forward", "adjoint"),
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

  If positive, print tracing information; higher values give more detail
  (see source for "L-BFGS-B" trace levels).

- factr:

  Convergence tolerance factor for "L-BFGS-B"; converges when the
  objective reduction is within this factor of machine tolerance
  (default 1e7, i.e. ~1e-8).

- pgtol:

  Tolerance on the projected gradient for "L-BFGS-B"; 0 (default)
  suppresses the check.

- abstol:

  Absolute x-value tolerance for "L-BFGS-B"; 0 (default) suppresses the
  check.

- reltol:

  Relative x-value tolerance for "L-BFGS-B"; 0 (default) suppresses the
  check.

- lmm:

  Number of BFGS updates retained in "L-BFGS-B" (default 5).

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

  Method for calculating the covariance. `"r"` (the default) uses
  nlmixr2's
  [`nlmixr2Hess()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2Hess.md)
  Hessian; `""` skips the covariance step.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits. One value drives, with a single
  consistent formula, the inner/outer optimizer convergence tolerance
  (`10^-sigdig`), the boundary check tolerance (`5*10^(-sigdig+1)`), and
  the ODE solver tolerances: the `rtol` exponent IS `sigdig` and `atol`
  sits three orders below, so `rtol = 10^-sigdig`,
  `atol = 10^(-sigdig-3)` for every solver (stiff, non-stiff or
  auto-switching). The sensitivity (`atolSens`/`rtolSens`) tolerances
  match the main solve (the outer gradient and covariance are built from
  them); the steady-state (`ssAtol`/`ssRtol`) tolerances run one order
  looser. Keying the optimizer to the same `10^-sigdig` means it
  converges to exactly the precision the solve supports. At the default
  `sigdig = 4` this is `atol = 1e-7`, `rtol = 1e-4`.

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
#>          OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -685.961 1157.916 1172.639       -575.958         3265.19        173.2675
#> 
#> ── Time (sec $time): ──
#> 
#>             setup optimize covariance preprocess postprocess table compress
#> elapsed 0.3732025 1.500666  6.229e-06      0.042       0.007 0.029    0.001
#>             other
#> elapsed 0.1101252
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI)
#> E0  -0.7441 0.2174 29.22 -0.7441 (-1.170, -0.3180)
#> Em    8.703  7.380 84.80     8.703 (-5.762, 23.17)
#> E50   4.311  2.684 62.26    4.311 (-0.9498, 9.572)
#> g     2.000  FIXED FIXED                     2.000
#>  
#>   Covariance Type ($covMethod): r
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0343     0 -0.389 -0.744
#> 2 1     0.0489     1 -1.13  -0.743
#> 3 1     0.0501     0 -0.389 -0.743
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$lbfgsb3c

fit2$lbfgsb3c
#> $par
#>         E0         Em        E50 
#> -0.7440591  8.7028559  4.3109578 
#> 
#> $grad
#> [1] -0.0013035589 -0.0026641811  0.0003459083
#> 
#> $value
#> [1] 575.958
#> 
#> $counts
#> [1] 12 12
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
#> 
#> $scaleC
#> [1] 0.002595759 0.031364338 0.028740636
#> 
#> $par.scaled
#>         E0         Em        E50 
#> -480.26594  260.53448   81.40733 
#> 
#> $hessian
#>               E0           Em          E50
#> E0   0.001349488  0.001140202 -0.003468655
#> Em   0.001140202  0.002448925 -0.006581618
#> E50 -0.003468655 -0.006581618  0.018346748
#> 
#> $cov.scaled
#>           E0        Em       E50
#> E0  7013.738  8311.419  4307.620
#> Em  8311.419 55371.787 21435.157
#> E50 4307.620 21435.157  8721.964
#> 
#> $r
#>                E0            Em          E50
#> E0   0.0006747439  0.0005701008 -0.001734327
#> Em   0.0005701008  0.0012244623 -0.003290809
#> E50 -0.0017343273 -0.0032908092  0.009173374
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
