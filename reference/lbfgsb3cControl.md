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

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

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
#> lPop -686.1126 1157.764 1172.488      -575.8822        13586.13        444.8447
#> 
#> ── Time (sec $time): ──
#> 
#>             setup optimize covariance preprocess postprocess table compress
#> elapsed 0.3837033 1.566815  8.582e-06      0.047       0.007 0.031    0.001
#>             other
#> elapsed 0.1014735
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI)
#> E0  -0.7174 0.2252 31.39 -0.7174 (-1.159, -0.2760)
#> Em    10.96  15.41 140.5     10.96 (-19.23, 41.16)
#> E50   5.065  4.781 94.38     5.065 (-4.305, 14.44)
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
#> 1 1     0.0343     0 -0.398 -0.717
#> 2 1     0.0489     1 -1.11  -0.716
#> 3 1     0.0501     0 -0.398 -0.716
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$lbfgsb3c

fit2$lbfgsb3c
#> $par
#>         E0         Em        E50 
#> -0.7173516 10.9632473  5.0650861 
#> 
#> $grad
#> [1]  1.906166e-07 -1.300116e-08 -1.253361e-06
#> 
#> $value
#> [1] 575.8822
#> 
#> $counts
#> [1] 23 23
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
#>        E0        Em       E50 
#> -469.9770  332.6033  107.6464 
#> 
#> $hessian
#>                E0            Em          E50
#> E0   0.0013513360  0.0008546165 -0.002924798
#> Em   0.0008546165  0.0014173389 -0.004414305
#> E50 -0.0029247983 -0.0044143052  0.013975726
#> 
#> $cov.scaled
#>            E0        Em      E50
#> E0   7524.644  22581.32  8707.16
#> Em  22581.315 241242.84 80923.55
#> E50  8707.160  80923.55 27668.54
#> 
#> $r
#>                E0            Em          E50
#> E0   0.0006756680  0.0004273082 -0.001462399
#> Em   0.0004273082  0.0007086695 -0.002207153
#> E50 -0.0014623992 -0.0022071526  0.006987863
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
