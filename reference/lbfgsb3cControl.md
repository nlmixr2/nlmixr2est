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
  stopFile = NULL,
  maxTime = NULL,
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
  sensMethod = c("default", "forward"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", ""),
  adjObf = TRUE,
  aboveStartTol = 0.01,
  ci = 0.95,
  sigdig = 3,
  sigdigTable = NULL,
  ...,
  condWarn = 1000
)
```

## Arguments

- trace:

  Ignored: the thread-safe C++ L-BFGS-B prints nothing; use \`print\`
  for the iteration log. Kept so existing calls still work.

- factr:

  Convergence tolerance factor for "L-BFGS-B"; converges when the
  objective reduction is within this factor of machine tolerance.
  \`NULL\` (default) uses \`10^(-sigdig-2) / .Machine\$double.eps\` (at
  least 10), two orders tighter than \`sigdig\`, as
  \`foceiControl(lbfgsFactr=)\` does.

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

  maximum number of objective function evaluations.

- returnLbfgsb3c:

  return the lbfgsb3c output instead of the nlmixr2 fit

- stopFile:

  Path of a file that, once created, stops the estimation gracefully:
  the run ends at its current state (for FOCEi, the lowest objective
  evaluated so far), removes the file, notes the stop in \`\$runInfo\`
  and returns the fit. \`NULL\` (default) does not poll.

- maxTime:

  Wall-clock limit in seconds, counted from the start of the estimation.
  Once it is reached the optimization ends as for \`stopFile\` and the
  fit is finalized from its current state; the covariance step and
  tables are not limited. \`NULL\` (default) sets no limit.

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
  \`rxControl()\`. The solver's adaptive steps set the objective's noise
  floor, which grows with \`rtol\` and the number of observations. For a
  non-stiff model, \`rxControl(method = "dop853")\` gives a much
  smoother objective than the default \`"liblsoda"\` at the same
  \`rtol\`.

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
  \`"fd"\` uses the legacy finite-difference behavior. Also gates the
  analytic moving-boundary correction for a modeled \`alag()\`/\`f()\`
  on a \`linCmt()\` compartment; set \`"fd"\` if that model infuses a
  dose into the lagged/scaled compartment (rxode2/rxode2#1236), or if
  the regimen also doses an \*unlagged/unscaled\* compartment alongside
  the lagged/scaled one – a common design for estimating \`f()\` from
  paired IV+oral data (rxode2/rxode2#1237).

- sensMethod:

  Method used to compute the ODE parameter sensitivities. \`"forward"\`
  uses the classic variational (forward) sensitivity ODEs; \`"default"\`
  is the same thing.

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

- aboveStartTol:

  Warn (in \`\$runInfo\`) when the final objective is above the
  objective at the initial estimates by more than this fraction of its
  magnitude (at least 0.1 objective units). Use \`NA\` or \`Inf\` to
  disable.

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits. One value drives, with a single
  consistent formula, the inner/outer optimizer convergence tolerance
  (`10^-sigdig`), the boundary check tolerance (`5*10^(-sigdig+1)`), and
  the ODE solver tolerances `rtol = max(10^(-2*sigdig), 1e-12)`,
  `atol = rtol/100` for every solver. An objective solved to relative
  error `rtol` locates its minimum only to about `sqrt(rtol)`, so this
  is the solve a `10^-sigdig` optimizer needs. The sensitivity
  (`atolSens`/`rtolSens`) tolerances match the main solve; the
  steady-state (`ssAtol`/`ssRtol`) tolerances run one order looser. At
  the default `sigdig = 3` this is `atol = 1e-8`, `rtol = 1e-6`.

- sigdigTable:

  Significant digits in the final output table. If not specified
  (\`NULL\`), it defaults to \`sigdig\`.

- ...:

  Ignored parameters

- condWarn:

  warn when the condition number of the covariance's correlation matrix
  exceeds this; `NA` or `Inf` turns the warning off.

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
#> → finding duplicate expressions in nlm pred-only...
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
#> lPop -658.2519 1185.625 1200.348      -589.8126        3935.758        207.5884
#> 
#> ── Time (sec $time): ──
#> 
#>             setup optimize covariance preprocess postprocess table compress
#> elapsed 0.1592299 1.480194 1.8581e-05       0.06       0.011 0.027    0.001
#>             other
#> elapsed 0.1405575
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>       Est.    SE  %RSE Back-transformed(95%CI)
#> E0  -0.520 0.110  21.1 -0.520 (-0.736, -0.305)
#> Em    7.95  3.96  49.7      7.95 (0.202, 15.7)
#> E50   4.38  1.63  37.1       4.38 (1.19, 7.58)
#> g     2.00 FIXED FIXED                    2.00
#>  
#>   Covariance Type ($covMethod): r
#>   Some strong fixed parameter correlations exist ($cor) :
#>      cor:Em,E0 cor:E50,E0 cor:E50,Em 
#>     0.486      0.597      0.979  
#>  
#> 
#>   Information about run found ($runInfo):
#>    • highly correlated estimates: Em/E50 (0.98) 
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0343     0 -0.467 -0.520
#> 2 1     0.0489     1 -0.986 -0.519
#> 3 1     0.0501     0 -0.467 -0.519
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$lbfgsb3c

fit2$lbfgsb3c
#> $par
#>         E0         Em        E50 
#> -0.5203094  7.9542417  4.3848589 
#> 
#> $grad
#>          E0          Em         E50 
#>  0.02293042  0.04393320 -0.02259589 
#> 
#> $value
#> [1] 589.8126
#> 
#> $counts
#> [1] 16 16
#> 
#> $convergence
#> [1] 0
#> 
#> $message
#> [1] "CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH"
#> 
#> $scaleC
#> [1] 2.0 2.0 0.5
#> 
#> $par.scaled
#>        E0        Em       E50 
#> -1.510155  2.727121  5.769718 
#> 
#> $hessian
#>            E0        Em        E50
#> E0  825.47069 55.892815 -41.565402
#> Em   55.89282 10.066091  -6.550727
#> E50 -41.56540 -6.550727   4.409652
#> 
#> $cov.scaled
#>              E0         Em        E50
#> E0  0.003018639 0.05279949  0.1068896
#> Em  0.052799488 3.91123200  6.3079914
#> E50 0.106889611 6.30799141 10.6051096
#> 
#> $r
#>            E0        Em        E50
#> E0  825.47069 55.892815 -41.565402
#> Em   55.89282 10.066091  -6.550727
#> E50 -41.56540 -6.550727   4.409652
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
