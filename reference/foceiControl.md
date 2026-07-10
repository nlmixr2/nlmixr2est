# Control Options for FOCEi

Control Options for FOCEi

## Usage

``` r
foceiControl(
  sigdig = 4,
  ...,
  epsilon = NULL,
  maxInnerIterations = 1000,
  maxOuterIterations = 5000,
  n1qn1nsim = NULL,
  print = 1L,
  printNcol = NULL,
  scaleTo = 1,
  scaleObjective = 0,
  normType = c("rescale2", "mean", "rescale", "std", "len", "constant"),
  scaleType = c("nlmixr2", "norm", "mult", "multAdd"),
  scaleCmax = 1e+05,
  scaleCmin = 1e-05,
  scaleC = NULL,
  scaleC0 = 1e+05,
  derivEps = rep(20 * sqrt(.Machine$double.eps), 2),
  derivMethod = c("switch", "forward", "central"),
  derivSwitchTol = NULL,
  covDerivMethod = c("central", "forward"),
  covMethod = c("analytic", "r,s", "r", "s", ""),
  covSolveTol = NULL,
  covFull = TRUE,
  fast = FALSE,
  hessEps = (.Machine$double.eps)^(1/3),
  hessEpsLlik = (.Machine$double.eps)^(1/3),
  optimHessType = c("central", "forward"),
  optimHessCovType = c("central", "forward"),
  censOption = c("gauss", "laplace"),
  eventType = c("central", "forward"),
  eventSens = c("jump", "fd"),
  centralDerivEps = rep(20 * sqrt(.Machine$double.eps), 2),
  lbfgsLmm = 7L,
  lbfgsPgtol = 0,
  lbfgsFactr = NULL,
  eigen = TRUE,
  diagXform = c("sqrt", "log", "identity"),
  iovXform = c("sd", "var", "logsd", "logvar"),
  sumProd = FALSE,
  optExpression = TRUE,
  literalFix = TRUE,
  literalFixRes = TRUE,
  ci = 0.95,
  useColor = NULL,
  boundTol = NULL,
  calcTables = TRUE,
  noAbort = TRUE,
  interaction = TRUE,
  foce = c("nonmem", "foce+"),
  cholSEtol = (.Machine$double.eps)^(1/3),
  cholAccept = 0.001,
  resetEtaP = 0.15,
  resetThetaP = 0.05,
  resetThetaFinalP = 0.15,
  diagOmegaBoundUpper = 5,
  diagOmegaBoundLower = 100,
  cholSEOpt = FALSE,
  cholSECov = FALSE,
  fo = FALSE,
  covTryHarder = FALSE,
  outerOpt = c("nlminb", "lbfgsb3c", "bobyqa", "L-BFGS-B", "mma", "lbfgsbLG", "slsqp",
    "uobyqa", "newuoa"),
  innerOpt = c("n1qn1", "BFGS"),
  rhobeg = 0.2,
  rhoend = NULL,
  npt = NULL,
  rel.tol = NULL,
  x.tol = NULL,
  eval.max = 4000,
  iter.max = 2000,
  abstol = NULL,
  reltol = NULL,
  resetHessianAndEta = FALSE,
  muModel = c("none", "irls", "lin"),
  muRefCovAlg = TRUE,
  muModelTol = 0.001,
  muModelMaxCycles = 10L,
  stateTrim = Inf,
  shi21maxOuter = 0L,
  shi21maxInner = 20L,
  shi21maxInnerCov = 20L,
  shi21maxFD = 20L,
  gillK = 10L,
  gillStep = 4,
  gillFtol = 0,
  gillRtol = sqrt(.Machine$double.eps),
  gillKcov = 10L,
  gillKcovLlik = 10L,
  gillStepCovLlik = 4.5,
  gillStepCov = 2,
  gillFtolCov = 0,
  gillFtolCovLlik = 0,
  rmatNorm = TRUE,
  rmatNormLlik = TRUE,
  smatNorm = TRUE,
  smatNormLlik = TRUE,
  covGillF = TRUE,
  optGillF = TRUE,
  covSmall = 1e-05,
  adjLik = TRUE,
  gradTrim = Inf,
  maxOdeRecalc = 5,
  odeRecalcFactor = 10^(0.5),
  gradCalcCentralSmall = 1e-04,
  gradCalcCentralLarge = 10000,
  etaNudge = qnorm(1 - 0.05/2)/sqrt(3),
  etaNudge2 = qnorm(1 - 0.05/2) * sqrt(3/5),
  nRetries = 3,
  seed = 42,
  resetThetaCheckPer = 0.1,
  etaMat = NULL,
  repeatGillMax = 1,
  stickyRecalcN = 4,
  indTolRelax = TRUE,
  gradProgressOfvTime = 10,
  addProp = c("combined2", "combined1"),
  badSolveObjfAdj = 100,
  compress = FALSE,
  rxControl = NULL,
  sigdigTable = NULL,
  fallbackFD = FALSE,
  smatPer = 0.6,
  sdLowerFact = 0.001,
  zeroGradFirstReset = TRUE,
  zeroGradRunReset = TRUE,
  zeroGradBobyqa = TRUE,
  mceta = -2L,
  warm = c("calc", "save"),
  nAGQ = 0,
  agqLow = -Inf,
  agqHi = Inf,
  sensMethod = c("default", "forward", "adjoint"),
  boundedTransform = TRUE
)
```

## Arguments

- sigdig:

  Optimization significant digits; controls the inner/outer optimization
  tolerance (`10^-sigdig`), ODE solver tolerance (`0.5*10^(-sigdig-2)`,
  or `0.5*10^(-sigdig-1.5)` for sensitivity/steady-state with liblsoda),
  and boundary check tolerance (`5*10^(-sigdig+1)`).

- ...:

  Ignored parameters

- epsilon:

  Precision of estimate for n1qn1 optimization.

- maxInnerIterations:

  Number of iterations for n1qn1 optimization.

- maxOuterIterations:

  Maximum number of L-BFGS-B optimization for outer problem.

- n1qn1nsim:

  Number of function evaluations for n1qn1 optimization.

- print:

  Either a scalar print-frequency (\`0\` = suppress, \`1\` (default) =
  every evaluation, \`N\` = every Nth), OR a pre-built
  \[iterPrintControl()\] object. Equivalent to \`iterPrintControl(every
  = print, ncol = printNcol, useColor = useColor)\`.

- printNcol:

  Integer (or \`NULL\`) parameter columns per row before wrapping.
  \`NULL\` (default) uses \`floor((getOption("width") - 23) / 12)\`.

- scaleTo:

  Scale the initial parameter estimate to this value. By default this
  is 1. When zero or below, no scaling is performed.

- scaleObjective:

  Scale the initial objective function to this value. By default this is
  0 (meaning do not scale)

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

- scaleC0:

  Number to adjust the scaling factor by if the initial gradient is
  zero.

- derivEps:

  Forward difference tolerances (relative, absolute); step size
  `h = abs(x)*derivEps[1] + derivEps[2]`.

- derivMethod:

  Derivative method for the outer problem: "switch", "central", or
  "forward". "switch" starts forward and toggles to central when
  `abs(delta(OFV)) <= derivSwitchTol`.

- derivSwitchTol:

  The tolerance to switch forward to central differences.

- covDerivMethod:

  indicates the method for calculating the derivatives while calculating
  the covariance components (Hessian and S).

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

- covSolveTol:

  absolute/relative ODE tolerance for the covariance solves – the
  augmented-sensitivity solves behind `covMethod="analytic"` and the
  perturbed solves behind the finite-difference methods. `NULL`
  (default) derives a tight tolerance from `sigdig`; supply a number to
  override it.

- covFull:

  controls the shape of `fit$cov`. `FALSE` (default) installs only the
  structural-theta block (the NONMEM-matched theta covariance, matching
  the historical finite-difference `fit$cov` shape for backwards
  compatibility); `TRUE` installs the full theta + residual sigma +
  Omega covariance – assembled analytically for `covMethod="analytic"`,
  or by central finite differences of the objective over the same
  parameter set for the finite-difference methods (perturbing Omega on
  the variance-covariance scale, with the per-parameter Gill (1983) step
  and the 5-point/4-point stencils that `foceiCalcR` uses). The theta
  standard errors are identical either way.

- fast:

  When `TRUE`, compute the outer (population) gradient analytically from
  Almquist (2015) sensitivity equations instead of by finite
  differences, and use the Eq-48 random-effect extrapolation for the
  next inner-problem starting values. Requires an analytic-scope model
  (single additive/proportional Gaussian endpoint); out-of-scope models
  fall back to the finite-difference gradient with a message (linCmt()
  and log-likelihood models downgrade to `fast=FALSE` up front). When
  unspecified, the outer optimizer defaults to `"lbfgsb3c"` (vs
  `"nlminb"` for `fast=FALSE`); pairing `fast=TRUE` with a
  derivative-free `outerOpt` reverts to `fast=FALSE`. The `*f` methods
  (e.g. `foceif`) default this to `TRUE`.

- hessEps:

  is a double value representing the epsilon for the Hessian
  calculation. This is used for the R matrix calculation.

- hessEpsLlik:

  is a double value representing the epsilon for the Hessian calculation
  when doing focei generalized log-likelihood estimation. This is used
  for the R matrix calculation.

- optimHessType:

  Hessian type for numeric-difference individual Hessians in generalized
  log-likelihood estimation: "central" (matches R's \`optimHess()\`,
  default) or "forward" (faster).

- optimHessCovType:

  Hessian type for numeric-difference individual Hessians used for the
  covariance step/final likelihood: "central" (more accurate, used here)
  or "forward".

- censOption:

  Treatment of the second derivative for censored (M2/M3/M4/BLQ)
  observations in the FOCEI family. `"gauss"` (the default) keeps the
  historic uncensored Gauss-Newton curvature, matching common PMx tools;
  `"laplace"` uses the exact censored second derivative of the objective
  (a proper Laplace inner Hessian and analytic covariance). Accepted by
  `saemControl`/`nlmControl` for a uniform interface but inert there –
  SAEM (stochastic EM) has no Laplace inner Hessian, and NLM uses a
  finite-difference Hessian that already reflects censoring exactly.

- eventType:

  Event gradient type for dosing events; Can be "central" or "forward"

- eventSens:

  Controls how dosing/event-parameter (\`alag\`, \`F\`, \`rate\`,
  \`dur\`) sensitivities are computed for THETA/ETA gradients:
  \`"jump"\` (default) uses rxode2's analytic event sensitivities;
  \`"fd"\` uses the legacy finite-difference behavior.

- centralDerivEps:

  Central difference tolerances (relative, absolute); step size
  `h = abs(x)*derivEps[1] + derivEps[2]`.

- lbfgsLmm:

  An integer giving the number of BFGS updates retained in the
  "L-BFGS-B" method, It defaults to 7.

- lbfgsPgtol:

  Projected-gradient convergence tolerance for "L-BFGS-B": iteration
  stops when `max(| proj g_i |) <= lbfgsPgtol`. Defaults to \`0\` (check
  suppressed).

- lbfgsFactr:

  Convergence factor for "L-BFGS-B": converges when the objective
  reduction is within `lbfgsFactr * .Machine$double.eps`. Default
  \`1e10\` (~4 sigdigs, `2e-6`).

- eigen:

  A boolean indicating if eigenvectors are calculated to include a
  condition number calculation.

- diagXform:

  Transformation used on the diagonal of `chol(solve(omega))` (the
  FOCEi-estimated parameters): one of `"sqrt"` (default), `"log"`, or
  `"identity"`.

- iovXform:

  Transformation used on the diagonal of the IOV: one of `"sd"`,
  `"var"`, `"logsd"`, or `"logvar"`.

- sumProd:

  Is a boolean indicating if the model should change multiplication to
  high precision multiplication and sums to high precision sums using
  the PreciseSums package. By default this is `FALSE`.

- optExpression:

  Optimize the rxode2 expression to speed up calculation. By default
  this is turned on.

- literalFix:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- literalFixRes:

  boolean, substitute fixed population values as literals and re-adjust
  ui and parameter estimates after optimization; Default is \`TRUE\`.

- ci:

  Confidence level for some tables. By default this is 0.95 or 95%
  confidence.

- useColor:

  Logical (or \`NULL\`) emit ANSI bold/color escapes in the iteration
  print. \`NULL\` (default) defers to \[crayon::has_color()\].

- boundTol:

  Tolerance for boundary issues.

- calcTables:

  This boolean is to determine if the foceiFit will calculate tables. By
  default this is `TRUE`

- noAbort:

  Boolean to indicate if you should abort the FOCEi evaluation if it
  runs into troubles. (default TRUE)

- interaction:

  Boolean indicate FOCEi should be used (TRUE) instead of FOCE (FALSE)

- foce:

  Controls how FOCE (`interaction = FALSE`) evaluates the residual
  variance R in the inner objective; ignored for FOCEi. Either
  `"nonmem"` (default) or `"foce+"`:

  - `"nonmem"` freezes R at the `eta = 0` population prediction and
    holds it constant across the inner optimization, matching NONMEM's
    FOCE. Advantage: reproduces NONMEM FOCE objective and standard
    errors, and an ODE model agrees with its closed-form (`linCmt`)
    equivalent. Disadvantage: R ignores the individual (conditional)
    heteroscedasticity, so it can be slightly less accurate than
    `"foce+"` for proportional/combined error.

  - `"foce+"` evaluates R at the current conditional `eta` (the live
    variance), keeping the truncated FOCE inner gradient. Advantage:
    uses the conditional variance and is a bit more accurate than
    NONMEM's FOCE in some cases. Disadvantage: does not match NONMEM
    FOCE. This was the FOCE behavior in nlmixr2est 6.0.1 and earlier.
    This does not use the gradient of `eta` like the full `focei`
    method, so it is not as accurate as `focei`.

- cholSEtol:

  tolerance for Generalized Cholesky Decomposition. Defaults to
  suggested (.Machine\$double.eps)^(1/3)

- cholAccept:

  Tolerance to accept a Generalized Cholesky Decomposition for a R or S
  matrix.

- resetEtaP:

  P-value for resetting an individual ETA to 0 during optimization,
  based on a z-test of `chol(omega^-1) %*% eta` or `eta/sd(allEtas)`.
  \`0\` = never reset, \`1\` = always reset.

- resetThetaP:

  P-value for resetting mu-referenced THETAs based on ETA drift, checked
  at the start and near a local minimum (see `resetThetaCheckPer`).
  \`0\` = never reset; \`1\` is not allowed.

- resetThetaFinalP:

  represents the p-value for reseting the population mu-referenced THETA
  parameters based on ETA drift during optimization, and resetting the
  optimization one final time.

- diagOmegaBoundUpper:

  Upper bound of the diagonal omega matrix, as
  `diag(omega)*diagOmegaBoundUpper`. \`1\` = no upper bound.

- diagOmegaBoundLower:

  Lower bound of the diagonal omega matrix, as
  `diag(omega)/diagOmegaBoundLower`. \`1\` = no lower bound.

- cholSEOpt:

  Boolean indicating if the generalized Cholesky should be used while
  optimizing.

- cholSECov:

  Boolean indicating if the generalized Cholesky should be used while
  calculating the Covariance Matrix.

- fo:

  is a boolean indicating if this is a FO approximation routine.

- covTryHarder:

  If the R matrix is non-positive definite and cannot be corrected to be
  non-positive definite try estimating the Hessian on the unscaled
  parameter space.

- outerOpt:

  optimization method for the outer problem

- innerOpt:

  optimization method for the inner problem (not implemented yet.)

- rhobeg:

  Initial trust region radius for the bobyqa outer optimizer (with
  \`rhoend\`, must satisfy \`0 \< rhoend \< rhobeg\`). Default \`0.2\`
  (20 \`abs(upper-lower)/2\`. (bobyqa)

- rhoend:

  Final trust region radius. If not defined, \`10^(-sigdig-1)\` is used.
  (bobyqa)

- npt:

  Number of points for bobyqa's quadratic approximation to the
  objective; must be in \`\[n+2, (n+1)(n+2)/2\]\`. Defaults to \`2\*n +
  1\`. (bobyqa)

- rel.tol:

  Relative tolerance before nlminb stops (nlmimb).

- x.tol:

  X tolerance for nlmixr2 optimizer

- eval.max:

  Number of maximum evaluations of the objective function (nlmimb)

- iter.max:

  Maximum number of iterations allowed (nlmimb)

- abstol:

  Absolute tolerance for nlmixr2 optimizer (BFGS)

- reltol:

  tolerance for nlmixr2 (BFGS)

- resetHessianAndEta:

  is a boolean representing if the individual Hessian is reset when ETAs
  are reset using the option `resetEtaP`.

- muModel:

  Selects the mu-referenced-FOCEI-family regression variant for
  theta/eta in a mu-ref covariate relationship (see `muRefCovAlg`):
  `"none"` (default, ordinary FOCEI); `"lin"`
  (`mufocei`/`mufoce`/`muagq`/ `mulaplace`: population theta and
  covariate coefficient(s) per mu-ref-covariate group are excluded from
  the outer optimizer and re-derived in C++ by closed-form OLS
  regression of each subject's back-calculated value on the
  covariate(s), residual becomes that subject's eta; repeats until
  convergence, see `muModelTol`/ `muModelMaxCycles`); or `"irls"`
  (`irlsfocei`/`irlsfoce`/`irlsagq`/`irlslaplace`: same mechanism,
  reweighted by inner-optimization curvature).

  A mu-ref-covariate theta with a finite bound falls back to ordinary
  bounded outer-optimizer handling with a warning (a bound on the
  group's population theta excludes the whole group; a bound on one
  covariate coefficient excludes only that covariate).

- muRefCovAlg:

  When \`TRUE\` (default), algebraic expressions that can be
  mu-referenced are internally rewritten as mu-referenced covariates and
  restored after optimization. Mirrors
  `saemControl(muRefCovAlg=)`/`nlmeControl(muRefCovAlg=)`; for
  `foceiControl()` only takes effect when `muModel != "none"`.

- muModelTol:

  Convergence tolerance for the mu-referenced-FOCEI-family "re-optimize
  etas, then regress" cycle (`muModel != "none"`): repeats until the max
  mu-group theta change drops below this value or `muModelMaxCycles` is
  reached.

- muModelMaxCycles:

  Maximum number of "re-optimize etas, regress" cycles per outer
  iteration (see `muModel`, `muModelTol`).

- stateTrim:

  Trim state amounts/concentrations to this value.

- shi21maxOuter:

  The maximum number of steps for the optimization of the
  forward-difference step size. When not zero, use this instead of Gill
  differences.

- shi21maxInner:

  The maximum number of steps for the optimization of the individual
  Hessian matrices in the generalized likelihood problem. When 0,
  un-optimized finite differences are used.

- shi21maxInnerCov:

  The maximum number of steps for the optimization of the individual
  Hessian matrices in the generalized likelihood problem for the
  covariance step. When 0, un-optimized finite differences are used.

- shi21maxFD:

  The maximum number of steps for the optimization of the forward
  difference step size when using dosing events (lag time, modeled
  duration/rate and bioavailability)

- gillK:

  Max steps to determine the optimal forward/central difference step
  size per parameter (Gill 1983). \`0\` = no optimal step size
  determined.

- gillStep:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration the new step size = (prior step size)\*gillStep

- gillFtol:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates.

- gillRtol:

  The relative tolerance used for Gill 1983 determination of optimal
  step size.

- gillKcov:

  Max steps to determine the optimal forward/central difference step
  size per parameter (Gill 1983) during the covariance step. \`0\` = no
  optimal step size determined.

- gillKcovLlik:

  Same as `gillK` but for the generalized focei log-likelihood method
  (Gill 1986).

- gillStepCovLlik:

  Same as above but during generalized focei log-likelihood

- gillStepCov:

  When looking for the optimal forward difference step size, this is
  This is the step size to increase the initial estimate by. So each
  iteration during the covariance step is equal to the new step size =
  (prior step size)\*gillStepCov

- gillFtolCov:

  The gillFtol is the gradient error tolerance that is acceptable before
  issuing a warning/error about the gradient estimates during the
  covariance step.

- gillFtolCovLlik:

  Same as above but applied during generalized log-likelihood
  estimation.

- rmatNorm:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the R matrix

- rmatNormLlik:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the R matrix if you are using generalized
  log-likelihood Hessian matrix.

- smatNorm:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the S matrix

- smatNormLlik:

  A parameter to normalize gradient step size by the parameter value
  during the calculation of the S matrix if you are using the
  generalized log-likelihood.

- covGillF:

  Use the Gill calculated optimal Forward difference step size for the
  instead of the central difference step size during the central
  difference gradient calculation.

- optGillF:

  Use the Gill calculated optimal Forward difference step size for the
  instead of the central difference step size during the central
  differences for optimization.

- covSmall:

  Small number used to compare covariance estimates (sandwich vs R/S
  matrix) before rejecting one as too small to be the final covariance
  estimate.

- adjLik:

  When \`TRUE\`, adjusts the likelihood by the 2\*pi constant nlmixr2's
  objective function otherwise omits (to match NONMEM), more closely
  matching nlme/SAS likelihood approximations. The objective function
  itself always matches NONMEM regardless.

- gradTrim:

  The parameter to adjust the gradient to if the \|gradient\| is very
  large.

- maxOdeRecalc:

  Maximum number of times to reduce the ODE tolerances and try to
  resolve the system if there was a bad ODE solve.

- odeRecalcFactor:

  The ODE recalculation factor when ODE solving goes bad, this is the
  factor the rtol/atol is reduced

- gradCalcCentralSmall:

  A small number that represents the value where \|grad\| \<
  gradCalcCentralSmall where forward differences switch to central
  differences.

- gradCalcCentralLarge:

  A large number that represents the value where \|grad\| \>
  gradCalcCentralLarge where forward differences switch to central
  differences.

- etaNudge:

  When n1qn1 optimization of an ETA (starting at zero) misbehaves, reset
  the Hessian and nudge the ETA up by this value, then down if it still
  doesn't move. Defaults to \`qnorm(1-0.05/2)\*1/sqrt(3)\`. Falls back
  to `etaNudge2`, then to zero (stop optimizing) if unsuccessful.

- etaNudge2:

  This is the second eta nudge. By default it is
  qnorm(1-0.05/2)\*sqrt(3/5), which is the n=3 quadrature point
  (excluding zero) times by the 0.95% normal region

- nRetries:

  If FOCEi doesn't fit with the current parameter estimates, randomly
  sample new parameter estimates and restart the problem. This is
  similar to 'PsN' resampling.

- seed:

  an object specifying if and how the random number generator should be
  initialized

- resetThetaCheckPer:

  represents objective function % percentage below which resetThetaP is
  checked.

- etaMat:

  Initial (or final) ETA estimates; can also be a prior fit, whose final
  ETAs are then used as initial values. By default, uses the last fit's
  ETAs if supplied, else all ETAs start at zero (\`NULL\`). \`NA\`
  disables reuse from a prior fit.

- repeatGillMax:

  If the tolerances were reduced when calculating the initial Gill
  differences, the Gill difference is repeated up to a maximum number of
  times defined by this parameter.

- stickyRecalcN:

  The number of bad ODE solves before reducing the atol/rtol for the
  rest of the problem.

- indTolRelax:

  When \`TRUE\` (default), only subjects whose ODE solve produced
  NaN/Inf have their tolerances relaxed, and the relaxed tolerance
  persists across optimizer calls (sticky). When \`FALSE\`, all subjects
  have their tolerances relaxed on each retry and tolerances are reset
  afterward.

- gradProgressOfvTime:

  This is the time for a single objective function evaluation (in
  seconds) to start progress bars on gradient evaluations

- addProp:

  Type of additive-plus-proportional error: \`"combined1"\`, where
  standard deviations add: \$\$y = f + (a + b\times f^c) \times
  \varepsilon\$\$; or \`"combined2"\`, where variances add: \$\$y = f +
  \sqrt{a^2 + b^2\times f^{2\times c}} \times \varepsilon\$\$. Here y =
  observed, f = predicted, a = additive sd, b = proportional/power sd, c
  = power exponent (1 in the proportional case).

- badSolveObjfAdj:

  The objective function adjustment when the ODE system cannot be
  solved. It is based on each individual bad solve.

- compress:

  Should the object have compressed items

- rxControl:

  \`rxode2\` ODE solving options during fitting, created with
  \`rxControl()\`

- sigdigTable:

  Significant digits in the final output table. If not specified, then
  it matches the significant digits in the \`sigdig\` optimization
  algorithm. If \`sigdig\` is NULL, use 3.

- fallbackFD:

  Fallback to the finite differences if the sensitivity equations do not
  solve.

- smatPer:

  Percentage of failed per-individual parameter gradients (replaced with
  the overall parameter gradient) out of the total (\`ntheta\*nsub\`)
  above which the S matrix is considered bad.

- sdLowerFact:

  Factor multiplying the estimate when the lower bound is zero for a
  standard-deviation error parameter (add.sd, prop.sd, etc); e.g.
  estimate 0.15 with lower bound 0 assumes a lower bound of 0.00015.
  \`0\` disables this.

- zeroGradFirstReset:

  When \`TRUE\` (default), reset a zero first gradient to
  \`sqrt(.Machine\$double.eps)\` instead of erroring; \`FALSE\` errors;
  \`NA\` ignores it only on the last reset attempt.

- zeroGradRunReset:

  When \`TRUE\` (default), reset a zero gradient encountered mid-run to
  \`sqrt(.Machine\$double.eps)\` instead of erroring.

- zeroGradBobyqa:

  When \`TRUE\` (default), a zero-gradient reset switches to the
  gradient-free bobyqa method; \`NA\` only does so for the first zero
  gradient.

- mceta:

  Monte Carlo sampling for the best initial ETA estimate (based on
  \`omega\`): \`-2\` (default) uses the Almquist (2015) Eq-48
  extrapolation \`eta^0 = eta\* + (d eta\*/d theta)(theta_new -
  theta_old)\` when the analytic gradient supplies \`d eta\*/d theta\`
  (\`fast = TRUE\`), accepting the extrapolated eta only when it is
  within the standardized-eta reset bound (else keeping the last eta, or
  resetting to 0 when that is also out of bound); \`-1\` jumps between
  the extrapolated eta and eta=0, keeping the better; both \`-2\` and
  \`-1\` fall back to keeping the last eta when no analytic \`d eta\*/d
  theta\` is available (\`fast = FALSE\`). \`0\` uses eta=0 for each
  inner optimization; for \`n\>0\`, the last eta, eta=0, and n-1 etas
  sampled from omega are each evaluated and the best (by inner
  objective) is used.

- warm:

  Seeding of the n1qn1 inner-optimization Hessian: \`"calc"\` (default)
  warm-starts each inner problem with the eta Hessian calculated for the
  objective function at the same starting eta, calculating it at the
  starting point when unavailable; \`"save"\` uses the classic
  self-initialized Hessian.

- nAGQ:

  Number of Gauss-Hermite adaptive quadrature points. \`0\` disables
  AGQ; \`1\` is equivalent to Laplace. Cost grows quickly with ETAs:
  once the EBE is found, expect \`nAGQ^neta\` (even \`nAGQ\`) or
  \`(nAGQ^neta)-1\` (odd \`nAGQ\`) additional evaluations per subject.

- agqLow:

  The lower bound for adaptive quadrature log-likelihood. By default
  this is -Inf; in the original nlmixr's gnlmm it was -700.

- agqHi:

  The upper bound for adaptive quadrature log-likelihood. By default
  this is Inf; in the original nlmixr's gnlmm was 400.

- sensMethod:

  Method used to compute the ODE parameter sensitivities: \`"default"\`
  (the default) defers to the global option
  \`getOption("nlmixr2est.adjoint")\` (itself \`"forward"\` by default);
  \`"forward"\` uses the classic variational (forward) sensitivity ODEs;
  \`"adjoint"\` uses the in-engine discrete adjoint with the matching
  adjoint (\`s\`) method.

- boundedTransform:

  When \`TRUE\` (default), bounded parameters are transformed for
  unbounded optimization methods and back-transformed for final
  estimates. \`FALSE\` optimizes on the original scale with bounds
  passed to the optimizer. \`NA\` transforms for optimization but skips
  the final back-transform.

## Value

The control object that changes the options for the FOCEi family of
estimation methods

## Details

Uses R's L-BFGS-B ([`optim`](https://rdrr.io/r/stats/optim.html)) for
the outer problem and BFGS
[`n1qn1`](https://nlmixr2.github.io/n1qn1c/reference/n1qn1.html)
(restoring the prior individual Hessian) for the inner problem, which is
left unscaled since eta estimates start near zero. The covariance step
is performed on the unscaled problem, so its condition number may differ
from the scaled problem's.

## References

Gill, P.E., Murray, W., Saunders, M.A., & Wright, M.H. (1983). Computing
Forward-Difference Intervals for Numerical Optimization. Siam Journal on
Scientific and Statistical Computing, 4, 310-321.

Shi, H.M., Xie, Y., Xuan, M.Q., & Nocedal, J. (2021). Adaptive
Finite-Difference Interval Estimation for Noisy Derivative-Free
Optimization.

## See also

[`optim`](https://rdrr.io/r/stats/optim.html)

[`n1qn1`](https://nlmixr2.github.io/n1qn1c/reference/n1qn1.html)

[`rxSolve`](https://nlmixr2.github.io/rxode2/reference/rxSolve.html)

Other Estimation control:
[`nlmixr2NlmeControl()`](https://nlmixr2.github.io/nlmixr2est/reference/nlmixr2NlmeControl.md),
[`saemControl()`](https://nlmixr2.github.io/nlmixr2est/reference/saemControl.md)

## Author

Matthew L. Fidler
