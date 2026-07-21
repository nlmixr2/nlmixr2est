# nlmixr2 defaults controls for nls

nlmixr2 defaults controls for nls

## Usage

``` r
nlsControl(
  maxiter = 10000,
  tol = NULL,
  minFactor = 1/1024,
  printEval = FALSE,
  warnOnly = FALSE,
  scaleOffset = 0,
  nDcentral = FALSE,
  algorithm = c("LM", "default", "plinear", "port"),
  ftol = NULL,
  ptol = NULL,
  gtol = 0,
  diag = list(),
  epsfcn = 0,
  factor = 100,
  maxfev = integer(),
  nprint = 0,
  solveType = c("grad", "fun"),
  stickyRecalcN = 4,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  indTolRelax = TRUE,
  eventType = c("central", "forward"),
  shiErr = (.Machine$double.eps)^(1/3),
  shi21maxFD = 20L,
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
  trace = FALSE,
  rxControl = NULL,
  optExpression = TRUE,
  sumProd = FALSE,
  literalFix = TRUE,
  returnNls = FALSE,
  addProp = c("combined2", "combined1"),
  eventSens = c("jump", "fd"),
  calcTables = TRUE,
  compress = TRUE,
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  ...
)
```

## Arguments

- maxiter:

  A positive integer specifying the maximum number of iterations
  allowed.

- tol:

  A positive numeric value specifying the tolerance level for the
  relative offset convergence criterion.

- minFactor:

  A positive numeric value specifying the minimum step-size factor
  allowed on any step in the iteration. The increment is calculated with
  a Gauss-Newton algorithm and successively halved until the residual
  sum of squares has been decreased or until the step-size factor has
  been reduced below this limit.

- printEval:

  a logical specifying whether the number of evaluations (steps in the
  gradient direction taken each iteration) is printed.

- warnOnly:

  a logical specifying whether
  [`nls()`](https://rdrr.io/r/stats/nls.html) should return instead of
  signalling an error in the case of termination before convergence.
  Termination before convergence happens upon completion of `maxiter`
  iterations, in the case of a singular gradient, and in the case that
  the step-size factor is reduced below `minFactor`.

- scaleOffset:

  a constant to be added to the denominator of the relative offset
  convergence criterion calculation to avoid a zero divide in the case
  where the fit of a model to data is very close. The default value of
  `0` keeps the legacy behaviour of
  [`nls()`](https://rdrr.io/r/stats/nls.html). A value such as `1` seems
  to work for problems of reasonable scale with very small residuals.

- nDcentral:

  only when *numerical* derivatives are used:
  [`logical`](https://rdrr.io/r/base/logical.html) indicating if
  *central* differences should be employed, i.e.,
  [`numericDeriv`](https://rdrr.io/r/stats/numericDeriv.html)`(*, central=TRUE)`
  be used.

- algorithm:

  character string specifying the algorithm to use. The default
  algorithm is a Gauss-Newton algorithm. Other possible values are
  `"plinear"` for the Golub-Pereyra algorithm for partially linear
  least-squares models and `"port"` for the ‘nl2sol’ algorithm from the
  Port library – see the references. Can be abbreviated.

- ftol:

  non-negative numeric. Termination occurs when both the actual and
  predicted relative reductions in the sum of squares are at most
  `ftol`. Therefore, `ftol` measures the relative error desired in the
  sum of squares.

- ptol:

  non-negative numeric. Termination occurs when the relative error
  between two consecutive iterates is at most `ptol`. Therefore, `ptol`
  measures the relative error desired in the approximate solution.

- gtol:

  non-negative numeric. Termination occurs when the cosine of the angle
  between result of `fn` evaluation \\fvec\\ and any column of the
  Jacobian is at most `gtol` in absolute value. Therefore, `gtol`
  measures the orthogonality desired between the function vector and the
  columns of the Jacobian.

- diag:

  a list or numeric vector containing positive entries that serve as
  multiplicative scale factors for the parameters. Length of `diag`
  should be equal to that of `par`. If not, user-provided `diag` is
  ignored and `diag` is internally set.

- epsfcn:

  (used if `jac` is not provided) is a numeric used in determining a
  suitable step for the forward-difference approximation. This
  approximation assumes that the relative errors in the functions are of
  the order of `epsfcn`. If `epsfcn` is less than the machine precision,
  it is assumed that the relative errors in the functions are of the
  order of the machine precision.

- factor:

  positive numeric, used in determining the initial step bound. This
  bound is set to the product of `factor` and the
  \\\|\code{diag}\*\code{par}\|\\ if nonzero, or else to `factor`
  itself. In most cases `factor` should lie in the interval (0.1,100).
  100 is a generally recommended value.

- maxfev:

  integer; termination occurs when the number of calls to `fn` has
  reached `maxfev`. Note that `nls.lm` sets the value of `maxfev` to
  `100*(length(par) + 1)` if `maxfev = integer()`, where `par` is the
  list or vector of parameters to be optimized.

- nprint:

  is an integer; set `nprint` to be positive to enable printing of
  iterates

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

- trace:

  logical value indicating if a trace of the iteration progress should
  be printed. Default is `FALSE`. If `TRUE` the residual (weighted)
  sum-of-squares, the convergence criterion and the parameter values are
  printed at the conclusion of each iteration. Note that
  [`format()`](https://rdrr.io/r/base/format.html) is used, so these
  mostly depend on
  [`getOption`](https://rdrr.io/r/base/options.html)`("digits")`. When
  the `"plinear"` algorithm is used, the conditional estimates of the
  linear parameters are printed after the nonlinear parameters. When the
  `"port"` algorithm is used the objective function value printed is
  half the residual (weighted) sum-of-squares.

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

- returnNls:

  logical; when TRUE, will return the nls object instead of the nlmixr
  object

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

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

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
  auto-switching). The sensitivity (`atolSens`/`rtolSens`) and
  steady-state (`ssAtol`/`ssRtol`) tolerances run one order looser.
  Keying the optimizer to the same `10^-sigdig` means it converges to
  exactly the precision the solve supports. At the default `sigdig = 4`
  this is `atol = 1e-7`, `rtol = 1e-4`.

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

  Additional optional arguments. None are used at present.

## Value

nls control object

## Author

Matthew L. Fidler

## Examples

``` r
# \donttest{

one.cmt <- function() {
  ini({
   tka <- 0.45
   tcl <- log(c(0, 2.7, 100))
   tv <- 3.45
   add.sd <- 0.7
 })
 model({
   ka <- exp(tka)
   cl <- exp(tcl)
   v <- exp(tv)
   linCmt() ~ add(add.sd)
 })
}

# Uses nlsLM from minpack.lm if available

fit1 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls", nlsControl(algorithm="LM"))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of nls model...
#> ✔ done
#> → calculate ∂(f)/∂(θ)
#> → finding duplicate expressions in nls gradient...
#> → optimizing duplicate expressions in nls gradient...
#> → finding duplicate expressions in nls pred-only...
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
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 2280

# Uses port and respect parameter boundaries
fit2 <- nlmixr(one.cmt, nlmixr2data::theo_sd, est="nls", nlsControl(algorithm="port"))
#>  
#>  
#>  
#>  
#> ℹ parameter labels from comments are typically ignored in non-interactive mode
#> ℹ Need to run with the source intact to parse comments
#> → loading into symengine environment...
#> → pruning branches (`if`/`else`) of nls model...
#> ✔ done
#> → calculate ∂(f)/∂(θ)
#> → finding duplicate expressions in nls gradient...
#> → optimizing duplicate expressions in nls gradient...
#> → finding duplicate expressions in nls pred-only...
#>  
#>  
#>  
#>  
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
#> → compress origData in nlmixr2 object, save 6584
#> → compress parHistData in nlmixr2 object, save 2576

# You can access the underlying nls object with `$nls`
fit2$nls
#> Nonlinear regression model
#>   model: 0 ~ nlmixr2est::.nlmixrNlsFunValGrad(DV, tka, rxBoundedTr.tcl,     tv)
#>    data: nlmixr2est::.nlmixrNlsData()
#>             tka rxBoundedTr.tcl              tv 
#>         -1.0097         -0.4351          1.0423 
#>  residual sum-of-squares: 249.7
#> 
#> Algorithm "port", convergence message: relative convergence (4)

# }
```
