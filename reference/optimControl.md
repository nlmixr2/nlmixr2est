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
  abstol = NULL,
  reltol = NULL,
  alpha = 1,
  beta = 0.5,
  gamma = 2,
  REPORT = NULL,
  warn.1d.NelderMead = TRUE,
  type = NULL,
  lmm = 5,
  factr = NULL,
  pgtol = 0,
  temp = 10,
  tmax = 10,
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1/3),
  shi21maxFD = 20L,
  solveType = c("grad", "fun"),
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
  returnOptim = FALSE,
  addProp = c("combined2", "combined1"),
  eventSens = c("jump", "fd"),
  sensMethod = c("default", "forward", "adjoint"),
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", "optim", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  boundedTransform = TRUE,
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

  helps control the convergence of the \`"L-BFGS-B"\` method. It is a
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

- solveType:

  controls whether \`optim\` uses nlmixr2's analytical gradients
  (event-related parameters like lag time/duration/rate/F use Shi2021
  finite differences instead). \`"gradient"\` supplies the gradient and
  lets \`optim\` compute the finite-difference Hessian; \`"fun"\` lets
  \`optim\` compute both by finite differences. Only applies to the
  gradient-based methods: "BFGS", "CG", "L-BFGS-B".

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

- returnOptim:

  logical; when TRUE this will return the optim list instead of the
  nlmixr2 fit object

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

- boundedTransform:

  When \`TRUE\` (default), bounded parameters are transformed for
  unbounded optimization methods and back-transformed for final
  estimates. \`FALSE\` optimizes on the original scale with bounds
  passed to the optimizer. \`NA\` transforms for optimization but skips
  the final back-transform.

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
#> lPop -684.2964 1159.581 1174.304      -576.7903        56372.43        479.4776
#> 
#> ── Time (sec value$time): ──
#> 
#>             setup  optimize covariance preprocess postprocess table compress
#> elapsed 0.3350544 0.5232404   6.58e-06      0.047       0.006 0.025    0.001
#>             other
#> elapsed 0.1156987
#> 
#> ── (value$parFixed or value$parFixedDf): ──
#> 
#>        Est.     SE  %RSE   Back-transformed(95%CI)
#> E0  -0.2947 0.2131 72.33 -0.2947 (-0.7124, 0.1231)
#> Em    17.22  31.65 183.8     17.22 (-44.81, 79.26)
#> E50   7.688  8.845 115.1     7.688 (-9.649, 25.02)
#> g     2.000  FIXED FIXED                     2.000
#>  
#>   Covariance Type (value$covMethod): r (optim)
#>   Censoring (value$censInformation): No censoring
#> 
#> ── Fit Data (object value is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0429     0 -0.557 -0.294
#> 2 1     0.0693     1 -0.850 -0.293
#> 3 1     0.0743     0 -0.557 -0.293
#> # ℹ 997 more rows
# }
```
