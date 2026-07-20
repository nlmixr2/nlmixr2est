# nlmixr2 defaults controls for nlm

nlmixr2 defaults controls for nlm

## Usage

``` r
nlmControl(
  typsize = NULL,
  fscale = 1,
  print.level = 0,
  ndigit = NULL,
  gradtol = 1e-06,
  stepmax = NULL,
  steptol = 1e-06,
  iterlim = 10000,
  check.analyticals = FALSE,
  returnNlm = FALSE,
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
  censOption = c("gauss", "laplace"),
  eventSens = c("jump", "fd"),
  sensMethod = c("default", "forward", "adjoint"),
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
  calcTables = TRUE,
  compress = FALSE,
  covMethod = c("r", "nlm", ""),
  adjObf = TRUE,
  ci = 0.95,
  sigdig = 4,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  ...
)
```

## Arguments

- typsize:

  an estimate of the size of each parameter at the minimum.

- fscale:

  an estimate of the size of `f` at the minimum.

- print.level:

  this argument determines the level of printing which is done during
  the minimization process. The default value of `0` means that no
  printing occurs, a value of `1` means that initial and final details
  are printed and a value of 2 means that full tracing information is
  printed.

- ndigit:

  the number of significant digits in the function `f`.

- gradtol:

  a positive scalar giving the tolerance at which the scaled gradient is
  considered close enough to zero to terminate the algorithm. The scaled
  gradient is a measure of the relative change in `f` in each direction
  `p[i]` divided by the relative change in `p[i]`.

- stepmax:

  a positive scalar which gives the maximum allowable scaled step
  length. `stepmax` is used to prevent steps which would cause the
  optimization function to overflow, to prevent the algorithm from
  leaving the area of interest in parameter space, or to detect
  divergence in the algorithm. `stepmax` would be chosen small enough to
  prevent the first two of these occurrences, but should be larger than
  any anticipated reasonable step.

- steptol:

  A positive scalar providing the minimum allowable relative step
  length.

- iterlim:

  a positive integer specifying the maximum number of iterations to be
  performed before the program is terminated.

- check.analyticals:

  a logical scalar specifying whether the analytic gradients and
  Hessians, if they are supplied, should be checked against numerical
  derivatives at the initial parameter values. This can help detect
  incorrectly formulated gradients or Hessians.

- returnNlm:

  is a logical that allows a return of the \`nlm\` object

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

- censOption:

  Treatment of the second derivative for censored (M2/M3/M4/BLQ)
  observations in the FOCEI family. `"gauss"` (the default) keeps the
  historic uncensored Gauss-Newton curvature, matching common PMx tools;
  `"laplace"` uses the exact censored second derivative of the objective
  (a proper Laplace inner Hessian and analytic covariance). Accepted by
  `saemControl`/`nlmControl` for a uniform interface but inert there –
  SAEM (stochastic EM) has no Laplace inner Hessian, and NLM uses a
  finite-difference Hessian that already reflects censoring exactly.

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

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- covMethod:

  "r" uses nlmixr2's \`nlmixr2Hess()\` for the hessian, or "nlm" uses
  the hessian from \`stats::nlm(.., hessian=TRUE)\`; defaults to "nlm"
  when using nlmixr2's hessian/gradient for solving.

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

  additional arguments to be passed to `f`.

## Value

nlm control object

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

fit2 <- nlmixr(mod, dsn, est="nlm")
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
#> ── nlmixr² log-likelihood nlm ──
#> 
#>           OBJF      AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -688.1882 1155.689 1170.412      -574.8444         2549012        156229.8
#> 
#> ── Time (sec $time): ──
#> 
#>             setup  optimize covariance preprocess postprocess table compress
#> elapsed 0.2456175 0.7658132  5.619e-06      0.037       0.004  0.02    0.008
#>             other
#> elapsed 0.0775637
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>        Est.    SE  %RSE Back-transformed(95%CI)
#> E0  -0.6131 4.608 751.6 -0.6131 (-9.645, 8.419)
#> Em    9.016 203.0  2251   9.016 (-388.8, 406.8)
#> E50   4.661 86.21  1849   4.661 (-164.3, 173.6)
#> g     2.000 FIXED FIXED                   2.000
#>  
#>   Covariance Type ($covMethod): r (nlm)
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     relative gradient is close to zero, current iterate is probably solution 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0358     0 -0.433 -0.613
#> 2 1     0.0521     0 -0.433 -0.612
#> 3 1     0.0610     0 -0.433 -0.612
#> # ℹ 997 more rows

# you can also get the nlm output with fit2$nlm

fit2$nlm
#> $minimum
#> [1] 574.8444
#> 
#> $estimate
#>         E0         Em        E50 
#> -0.6130711  9.0163993  4.6613852 
#> 
#> $gradient
#> [1] -5.406938e-07 -1.822863e-06  3.341858e-06
#> 
#> $hessian
#>               E0           Em          E50
#> E0   0.001540707  0.001102975 -0.003277943
#> Em   0.001102975  0.002027278 -0.005561242
#> E50 -0.003277943 -0.005561242  0.013715325
#> 
#> $code
#> [1] 1
#> 
#> $iterations
#> [1] 17
#> 
#> $scaleC
#> [1] 0.002771929 0.032196427 0.029707265
#> 
#> $estimate.scaled
#>         E0         Em        E50 
#> -402.55106  263.51380   90.58702 
#> 
#> $cov.scaled
#>           E0       Em      E50
#> E0   2763643 10470185  4820804
#> Em  10470185 39741386 18293564
#> E50  4820804 18293564  8421378
#> 
#> $r
#>                E0            Em          E50
#> E0   0.0007703535  0.0005514877 -0.001638971
#> Em   0.0005514877  0.0010136390 -0.002780621
#> E50 -0.0016389714 -0.0027806211  0.006857663
#> 

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
