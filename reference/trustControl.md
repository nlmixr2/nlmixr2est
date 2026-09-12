# Control for the trust estimation method in nlmixr2

\`est="trust"\` is a trust-region Newton method for the population theta
vector, backed by the \`RcppTrust\` package's thread-safe
\`trust_solve_c()\`. Unlike every other nlm-family method
(\`nlm\`/\`nlminb\`/\`bobyqa\`/\`newuoa\`/
\`uobyqa\`/\`n1qn1\`/\`lbfgsb3c\`/\`optim\`), whose optimization loop
lives in R (one \`.Call()\` per iteration), \`trust\`'s entire loop runs
inside a single \`.Call()\` – \`RcppTrust\` needs no R API, so there is
no per-iteration R round-trip. It also supplies a full
analytic-or-Shi21-finite-difference gradient AND a full
finite-difference-of-the-gradient Hessian every iteration (there is no
analytic outer-theta Hessian in this package), so each outer iteration
costs roughly \`ntheta\` extra full population-gradient solves on top of
the gradient solve itself – state this cost plainly, it is the price of
true Newton-trust behavior rather than a derivative-free or quasi-Newton
method. \`trust\` is unbounded, like \`n1qn1\`/\`nlm\`: box constraints
are handled upstream by \`preProcessBoundedTransform.R\`, not by the
optimizer itself.

## Usage

``` r
trustControl(
  rinit = NULL,
  rmax = NULL,
  iterlim = 1000L,
  fterm = NULL,
  mterm = NULL,
  optimHessType = 1L,
  shi21maxHess = 20L,
  hessErr = NULL,
  hessianMethod = c("sr1", "fd", "bfgs", "bofill"),
  returnTrust = FALSE,
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
  sigdig = 3,
  sigdigTable = NULL,
  boundedTransform = TRUE,
  ...
)
```

## Arguments

- rinit:

  Initial trust-region radius, in the SAME scaled-parameter space every
  nlm-family method (bobyqa included) optimizes in. \`NULL\` (default)
  derives it from the scaled starting vector using \`minqa::bobyqa()\`'s
  own default-\`rhobeg\` formula, \`min(0.95, 0.2\*max(abs(par.ini)))\`,
  so swapping \`est="bobyqa"\` for \`est="trust"\` on the same model
  starts from a comparable trust-region size.

- rmax:

  Maximum trust-region radius, same scaled-parameter space. \`NULL\`
  (default) derives it as \`8 \* rinit\` – the same growth-ceiling
  multiplier already used for \`foceiControl(innerOpt="trust")\`'s per-
  subject eta trust region.

- iterlim:

  Maximum number of \`trust_solve_c()\` iterations.

- fterm, mterm:

  \`trust_solve_c()\`'s function-value and predicted- decrease
  convergence tolerances. \`NULL\` (default) uses \`10^(-sigdig-2)\`,
  two orders tighter than \`bobyqaControl()\`'s \`rhoend\` and
  \`foceiControl()\`'s inner \`epsilon\` (both the plain
  \`10^(-sigdig)\`), matching \`foceiControl(trustFterm=, trustMterm=)\`
  instead – the analogous tolerance for the other \`RcppTrust\`-backed
  solve in this package. \`mterm\` defaults to \`fterm\` when not given
  separately.

- optimHessType:

  Finite-difference type for the per-iteration outer Hessian: \`1L\`
  forward (default), \`2L\` central. Exposed explicitly (unlike
  bobyqa/n1qn1, where this is only ever set implicitly) because
  \`trust\` recomputes this Hessian every outer iteration, not once
  post-fit, making its accuracy/cost tradeoff far more consequential; a
  stiff or noisy model may benefit from \`optimHessType=2\`.

- shi21maxHess:

  Maximum Shi (2021) adaptive-step-size iterations for the per-iteration
  outer Hessian.

- hessErr:

  Target relative error for the per-iteration outer Hessian's Shi (2021)
  step-size search. \`NULL\` (default) uses
  \`(.Machine\$double.eps)^(1/3)\`, the same fallback \`.nlmSetupEnv()\`
  itself applies for every nlm-family method.

- hessianMethod:

  How the per-iteration outer Hessian is built. \`"sr1"\` (default) is
  the Symmetric Rank-1 quasi-Newton update (Nocedal & Wright,
  \*Numerical Optimization\*, 2nd ed., 2006, Eq. 6.24; Murtagh &
  Sargent, \*Comput. J.\* 13, 1970), built from consecutive outer
  iterations' gradients (already computed regardless of
  \`hessianMethod\`, so this adds no extra evaluations) – Nocedal &
  Wright's own recommendation for trust-region methods specifically,
  since (unlike BFGS) it is not forced positive definite, so it can
  represent indefinite curvature. \`"fd"\` instead recomputes the
  Hessian from scratch every outer iteration via \`nlmCalcHessian()\`'s
  Shi (2021) finite-difference-of-the-gradient (the
  \`optimHessType\`/\`shi21maxHess\`/ \`hessErr\` parameters above only
  apply to this method). \`"bfgs"\` is the damped BFGS update (Nocedal &
  Wright Procedure 18.2), always positive definite. \`"bofill"\` is
  Bofill's SR1/PSB blend (\*J. Comput. Chem.\* 15, 1-11, 1994). Every
  method seeds from one \`"fd"\`-style Hessian on the first outer
  iteration.

  \`"sr1"\` was originally made default from a benchmark
  (\`inst/benchmarks/benchmark-trust-outer.R\`) showing it ran faster
  with the same or slightly better accuracy than \`"fd"\`; that
  benchmark predated fixes for two real correctness bugs (issues \#994
  and \#996) that independently distorted several of its models' results
  for EVERY \`hessianMethod\` value alike, both upstream of Hessian
  construction (a \`scaleC\` blowup for a near-zero starting gradient;
  \`est="trust"\` missing the \`linCmt()\`-to-ODE translation another
  nlm-family method needs) – so it was briefly reverted to \`"fd"\`
  pending confirmation. Re-run after both fixes, \`"sr1"\`/\`"bofill"\`
  track \`"fd"\` closely (median \|objective diff\| vs \`bobyqa\`
  1.53/1.55 vs \`"fd"\`'s 1.55 across the corpus) and \`"bfgs"\` if
  anything tracks it slightly better (0.43) – confirming the earlier
  small \`"sr1"\`-vs-\`"fd"\` accuracy gap was, at least in part, noise
  from those two bugs, not a genuine difference between Hessian
  constructions for this OUTER problem, so \`"sr1"\` is restored as the
  default (faster, with no demonstrated accuracy cost for this problem).

  This OUTER problem does not have the failure mode that keeps the
  analogous inner-problem option (\`foceiControl(hessianMethod=)\`)
  defaulted to \`"fd"\`: that inner Hessian's log-determinant feeds
  directly into the reported per-subject objective (not just the step),
  and a quasi-Newton estimate was shown to bias it on a real PK model.
  \`nlmTrustObjfun()\`'s reported value (\`src/nlm.cpp\`) here is
  instead the plain log-likelihood, set before the Hessian is even
  touched, so a less-accurate \`hessianMethod\` can only cost step
  quality/convergence speed, not silently bias the reported number –
  which is what the re-benchmark above confirms in practice.

- returnTrust:

  return the raw \`nlmTrustFit()\` output list instead of the nlmixr2
  fit.

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
  \`"fd"\` uses the legacy finite-difference behavior. Also gates the
  analytic moving-boundary correction for a modeled \`alag()\`/\`f()\`
  on a \`linCmt()\` compartment; set \`"fd"\` if that model infuses a
  dose into the lagged/scaled compartment (rxode2/rxode2#1236), or if
  the regimen also doses an \*unlagged/unscaled\* compartment alongside
  the lagged/scaled one – a common design for estimating \`f()\` from
  paired IV+oral data (rxode2/rxode2#1237).

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- compress:

  Should the object have compressed items

- covMethod:

  Method for calculating the covariance. \`"r"\` (the default) reuses
  the LAST outer iteration's already-computed Hessian (skipping
  \`nlmixr2est\`'s own post-fit finite-difference Hessian recompute,
  since \`trust\` already has one in hand); \`""\` skips the covariance
  step.

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
  `sigdig = 3` this is `atol = 1e-6`, `rtol = 1e-3`.

- sigdigTable:

  Significant digits in the final output table. If not specified
  (\`NULL\`), it defaults to \`sigdig\`.

- boundedTransform:

  When \`TRUE\` (default), bounded parameters are transformed for
  unbounded optimization methods and back-transformed for final
  estimates. \`FALSE\` optimizes on the original scale with bounds
  passed to the optimizer. \`NA\` transforms for optimization but skips
  the final back-transform.

- ...:

  Ignored parameters

## Value

trust control structure

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

fit <- nlmixr(mod, dsn, est="trust")
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

print(fit)
#> ── nlmixr² log-likelihood trust ──
#> 
#>           OBJF     AIC      BIC Log-likelihood Condition#(Cov) Condition#(Cor)
#> lPop -679.7168 1164.16 1178.884      -579.0801        380.6064        65.74402
#> 
#> ── Time (sec $time): ──
#> 
#>             setup  optimize covariance preprocess postprocess table compress
#> elapsed 0.1326775 0.1715734  5.328e-06      0.042       0.006 0.017    0.001
#>              other
#> elapsed 0.07274377
#> 
#> ── ($parFixed or $parFixedDf): ──
#> 
#>       Est.    SE  %RSE Back-transformed(95%CI)
#> E0  -0.806 0.126  15.6  -0.806 (-1.05, -0.560)
#> Em    5.49  1.18  21.6       5.49 (3.17, 7.81)
#> E50   2.90 0.622  21.5       2.90 (1.68, 4.12)
#> g     2.00 FIXED FIXED                    2.00
#>  
#>   Covariance Type ($covMethod): r
#>   Some strong fixed parameter correlations exist ($cor) :
#>      cor:Em,E0 cor:E50,E0 cor:E50,Em 
#>     0.360      0.620      0.917  
#>  
#> 
#>   Censoring ($censInformation): No censoring
#>   Minimization message ($message):  
#>     converged 
#> 
#> ── Fit Data (object is a modified tibble): ──
#> # A tibble: 1,000 × 5
#>   ID      TIME    DV  IPRED      v
#>   <fct>  <dbl> <dbl>  <dbl>  <dbl>
#> 1 1     0.0383     0 -0.369 -0.805
#> 2 1     0.0509     0 -0.370 -0.805
#> 3 1     0.0543     1 -1.17  -0.804
#> # ℹ 997 more rows

# you can also get the raw trust output with fit$trust

fit$trust
#> $par
#>         E0         Em        E50 
#> -0.8062205  5.4865859  2.8981420 
#> 
#> $fval
#> [1] 579.0801
#> 
#> $hessian
#>            E0        Em       E50
#> E0  745.41393 104.57468 -68.97597
#> Em  104.57468  32.68177 -17.53328
#> E50 -68.97597 -17.53328  10.45523
#> 
#> $convergence
#> [1] 0
#> 
#> $iterations
#> [1] 18
#> 
#> $message
#> [1] "converged"
#> 
#> $scaleC
#> [1] 2.0 2.0 0.5
#> 
#> $par.scaled
#>        E0        Em       E50 
#> -1.653110  1.493293  2.796284 
#> 
#> $cov.scaled
#>              E0        Em        E50
#> E0  0.003955674 0.0133884 0.04854882
#> Em  0.013388401 0.3503137 0.67579822
#> E50 0.048548822 0.6757982 1.54924023
#> 
#> $r
#>            E0        Em       E50
#> E0  745.41393 104.57468 -68.97597
#> Em  104.57468  32.68177 -17.53328
#> E50 -68.97597 -17.53328  10.45523
#> 
# }
```
