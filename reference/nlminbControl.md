# nlmixr2 nlminb defaults

nlmixr2 nlminb defaults

## Usage

``` r
nlminbControl(
  eval.max = 200,
  iter.max = 150,
  trace = 0,
  abs.tol = 0,
  rel.tol = NULL,
  x.tol = NULL,
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
  addProp = c("combined2", "combined1"),
  eventSens = c("jump", "fd"),
  sensMethod = c("default", "forward", "adjoint"),
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

  Relative tolerance. When \`NULL\` (default) it is derived from
  \`sigdig\` the way \`foceiControl()\` does (\`10^(-sigdig-1)\`).

- x.tol:

  X tolerance. When \`NULL\` (default) it is derived from \`sigdig\`
  (\`10^(-sigdig-1)\`).

- xf.tol:

  false convergence tolerance. Defaults to \`2.2e-14\`.

- step.min:

  Minimum step size. Default to \`1.\`.

- step.max:

  Maximum step size. Default to \`1.\`.

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

  controls whether \`nlm\` uses nlmixr2's analytical gradients
  (event-related parameters like lag time/duration/rate/F use Shi2021
  finite differences instead): \`"hessian"\` builds a Hessian from the
  analytical gradient via finite differences, \`"gradient"\` supplies
  the gradient and lets \`nlm\` compute the finite-difference Hessian,
  and \`"fun"\` lets \`nlm\` compute both by finite differences.

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

  Hessian type for numeric-difference individual Hessians in generalized
  log-likelihood estimation: "central" (matches R's \`optimHess()\`,
  default) or "forward" (faster).

- hessErr:

  This represents the epsilon when optimizing the Hessian step size
  using the Shi2021 method.

- shi21maxHess:

  Maximum number of times to optimize the best step size for the hessian
  calculation

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
  Hessian; `"nlminb"` uses the optimizer's own Hessian; `""` skips the
  covariance step.

- adjObf:

  is a boolean to indicate if the objective function should be adjusted
  to be closer to NONMEM's default objective function. By default this
  is `TRUE`

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), the boundary check tolerance
  (`5*10^(-sigdig+1)`), and the ODE solver tolerances. The solver
  tolerances are split by solver stiffness and keep `atol` well below
  `rtol`: a stiff solver (`lsoda`/`liblsoda` – and any auto-switching
  method – plus `indLin`, the default) uses `rtol = 10^(-sigdig-3)`,
  `atol = 10^(-sigdig-5)`, while the non-stiff explicit `dop853` uses
  the looser `rtol = 10^-sigdig`, `atol = 10^(-sigdig-3)`. The
  sensitivity (`atolSens`/`rtolSens`) and steady-state
  (`ssAtol`/`ssRtol`) tolerances run one order looser than the
  corresponding main tolerance. At the default `sigdig = 4` a stiff
  solve is `atol = 1e-9`, `rtol = 1e-7`.

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
#> → compress origData in nlmixr2 object, save 8328
#> → compress parHistData in nlmixr2 object, save 2608

print(fit2)
#> ── nlmixr² log-likelihood nlminb ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -717.7048 1126.172 1140.896      -560.0862        143.6999        39.92278
#> 
#> ── Time (sec $time): ──
#> 
#>             setup  optimize covariance preprocess postprocess table compress
#> elapsed 0.3409622 0.6447436  7.343e-06      0.056       0.006  0.38    0.012
#>             other
#> elapsed 0.1102869
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>       Est.     SE  %RSE  Back-transformed(95%CI)
#> E0  -1.047 0.2689 25.67 -1.047 (-1.574, -0.5203)
#> Em   4.917  1.343 27.30     4.917 (2.286, 7.549)
#> E50  2.209 0.7169 32.46    2.209 (0.8038, 3.614)
#> g    2.000  FIXED FIXED                    2.000
#>  
#>   Covariance Type ($covMethod): r (nlminb)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     relative convergence (4) 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED     v
#>   <fct>  <dbl> <dbl>  <dbl> <dbl>
#> 1 1     0.0175     0 -0.301 -1.05
#> 2 1     0.0480     1 -1.35  -1.05
#> 3 1     0.0529     0 -0.301 -1.04
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$nlminb

fit2$nlminb
#> $par
#>        E0        Em       E50 
#> -1.047335  4.917202  2.208856 
#> 
#> $objective
#> [1] 560.0862
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 7
#> 
#> $evaluations
#> function gradient 
#>       14        8 
#> 
#> $message
#> [1] "relative convergence (4)"
#> 
#> $scaleC
#> [1] 0.002625647 0.033446592 0.031305839
#> 
#> $par.scaled
#>          E0          Em         E50 
#> -590.315830  131.067335    7.671457 
#> 
#> $hessian
#>               E0           Em         E50
#> E0   0.001316937  0.003308724 -0.00851920
#> Em   0.003308724  0.016106760 -0.03268196
#> E50 -0.008519200 -0.032681956  0.07905771
#> 
#> $cov.scaled
#>             E0        Em       E50
#> E0  10486.8694  860.6009 1485.8243
#> Em    860.6009 1611.3108  758.8434
#> E50  1485.8243  758.8434  524.4083
#> 
#> $r
#>                E0           Em         E50
#> E0   0.0006584684  0.001654362 -0.00425960
#> Em   0.0016543618  0.008053380 -0.01634098
#> E50 -0.0042596001 -0.016340978  0.03952885
#> 
# }
```
