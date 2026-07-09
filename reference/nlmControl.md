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
#> Error in (function (typsize = NULL, fscale = 1, print.level = 0, ndigit = NULL,     gradtol = 1e-06, stepmax = NULL, steptol = 1e-06, iterlim = 10000,     check.analyticals = FALSE, returnNlm = FALSE, solveType = c("hessian",         "grad", "fun"), stickyRecalcN = 4, maxOdeRecalc = 5,     odeRecalcFactor = 10^(0.5), indTolRelax = TRUE, eventType = c("central",         "forward"), shiErr = (.Machine$double.eps)^(1/3), shi21maxFD = 20L,     optimHessType = c("central", "forward"), hessErr = (.Machine$double.eps)^(1/3),     shi21maxHess = 20L, censOption = c("gauss", "laplace"), eventSens = c("jump",         "fd"), sensMethod = c("default", "forward", "adjoint"),     useColor = NULL, printNcol = NULL, print = 1L, normType = c("rescale2",         "mean", "rescale", "std", "len", "constant"), scaleType = c("nlmixr2",         "norm", "mult", "multAdd"), scaleCmax = 1e+05, scaleCmin = 1e-05,     scaleC = NULL, scaleTo = 1, gradTo = 1, rxControl = NULL,     optExpression = TRUE, sumProd = FALSE, literalFix = TRUE,     literalFixRes = TRUE, addProp = c("combined2", "combined1"),     calcTables = TRUE, compress = FALSE, covMethod = c("r", "nlm",         ""), adjObf = TRUE, ci = 0.95, sigdig = 4, sigdigTable = NULL,     boundedTransform = TRUE, ...) {    checkmate::assertNumeric(shiErr, lower = 0, any.missing = FALSE,         len = 1)    checkmate::assertNumeric(hessErr, lower = 0, any.missing = FALSE,         len = 1)    checkmate::assertIntegerish(shi21maxFD, lower = 1, any.missing = FALSE,         len = 1)    checkmate::assertIntegerish(shi21maxHess, lower = 1, any.missing = FALSE,         len = 1)    checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)    checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)    checkmate::assertLogical(literalFixRes, len = 1, any.missing = FALSE)    checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)    checkmate::assertNumeric(stepmax, lower = 0, len = 1, null.ok = TRUE,         any.missing = FALSE)    checkmate::assertIntegerish(print.level, lower = 0, upper = 2,         any.missing = FALSE)    checkmate::assertNumeric(ndigit, lower = 0, len = 1, any.missing = FALSE,         null.ok = TRUE)    checkmate::assertNumeric(gradtol, lower = 0, len = 1, any.missing = FALSE)    checkmate::assertNumeric(steptol, lower = 0, len = 1, any.missing = FALSE)    checkmate::assertIntegerish(iterlim, lower = 1, len = 1,         any.missing = FALSE)    checkmate::assertLogical(check.analyticals, len = 1, any.missing = FALSE)    checkmate::assertLogical(returnNlm, len = 1, any.missing = FALSE)    checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)    checkmate::assertLogical(compress, len = 1, any.missing = TRUE)    checkmate::assertLogical(adjObf, len = 1, any.missing = TRUE)    checkmate::assertLogical(boundedTransform, len = 1, any.missing = FALSE)    .xtra <- list(...)    .bad <- names(.xtra)    .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]    if (length(.bad) > 0) {        stop("unused argument: ", paste(paste0("'", .bad, "'",             sep = ""), collapse = ", "), call. = FALSE)    }    checkmate::assertIntegerish(stickyRecalcN, any.missing = FALSE,         lower = 0, len = 1)    checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE,         len = 1)    checkmate::assertNumeric(odeRecalcFactor, len = 1, lower = 1,         any.missing = FALSE)    checkmate::assertLogical(indTolRelax, any.missing = FALSE,         len = 1)    .genRxControl <- FALSE    if (!is.null(.xtra$genRxControl)) {        .genRxControl <- .xtra$genRxControl    }    if (is.null(ndigit)) {        ndigit <- sigdig    }    if (is.null(rxControl)) {        if (!is.null(sigdig)) {            rxControl <- rxode2::rxControl(sigdig = sigdig)        }        else {            rxControl <- rxode2::rxControl(atol = 1e-04, rtol = 1e-04)        }        .genRxControl <- TRUE    }    else if (inherits(rxControl, "rxControl")) {    }    else if (is.list(rxControl)) {        rxControl <- do.call(rxode2::rxControl, rxControl)    }    else {        stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'",             call = FALSE)    }    if (!is.null(sigdig)) {        checkmate::assertNumeric(sigdig, lower = 1, finite = TRUE,             any.missing = TRUE, len = 1)        if (is.null(sigdigTable)) {            sigdigTable <- round(sigdig)        }    }    if (is.null(sigdigTable)) {        sigdigTable <- 3    }    checkmate::assertIntegerish(sigdigTable, lower = 1, len = 1,         any.missing = FALSE)    .solveTypeIdx <- c(hessian = 3L, grad = 2L, fun = 1L)    if (checkmate::testIntegerish(solveType, len = 1, lower = 1,         upper = 6, any.missing = FALSE)) {        solveType <- as.integer(solveType)    }    else {        solveType <- setNames(.solveTypeIdx[match.arg(solveType)],             NULL)    }    if (missing(covMethod) && any(solveType == 2:3)) {        covMethod <- "nlm"    }    else {        covMethod <- match.arg(covMethod)    }    .eventTypeIdx <- c(central = 2L, forward = 1L)    if (checkmate::testIntegerish(eventType, len = 1, lower = 1,         upper = 6, any.missing = FALSE)) {        eventType <- as.integer(eventType)    }    else {        eventType <- setNames(.eventTypeIdx[match.arg(eventType)],             NULL)    }    .optimHessTypeIdx <- c(central = 2L, forward = 1L)    if (checkmate::testIntegerish(optimHessType, len = 1, lower = 1,         upper = 6, any.missing = FALSE)) {        optimHessType <- as.integer(optimHessType)    }    else {        optimHessType <- setNames(.optimHessTypeIdx[match.arg(optimHessType)],             NULL)    }    if (checkmate::testIntegerish(censOption, len = 1, lower = 0,         upper = 1, any.missing = FALSE)) {        censOption <- as.integer(censOption)    }    else {        censOption <- setNames(c(gauss = 0L, laplace = 1L)[match.arg(censOption)],             NULL)    }    eventSens <- match.arg(eventSens)    sensMethod <- match.arg(sensMethod)    eventSens <- match.arg(eventSens)    sensMethod <- match.arg(sensMethod)    .iterPrintControl <- .absorbIterPrintControl(print = print,         printNcol = printNcol, useColor = useColor, iterPrintControl = .xtra$iterPrintControl)    if (checkmate::testIntegerish(scaleType, len = 1, lower = 1,         upper = 4, any.missing = FALSE)) {        scaleType <- as.integer(scaleType)    }    else {        .scaleTypeIdx <- c(norm = 1L, nlmixr2 = 2L, mult = 3L,             multAdd = 4L)        scaleType <- setNames(.scaleTypeIdx[match.arg(scaleType)],             NULL)    }    .normTypeIdx <- c(rescale2 = 1L, rescale = 2L, mean = 3L,         std = 4L, len = 5L, constant = 6L)    if (checkmate::testIntegerish(normType, len = 1, lower = 1,         upper = 6, any.missing = FALSE)) {        normType <- as.integer(normType)    }    else {        normType <- setNames(.normTypeIdx[match.arg(normType)],             NULL)    }    checkmate::assertNumeric(scaleCmax, lower = 0, any.missing = FALSE,         len = 1)    checkmate::assertNumeric(scaleCmin, lower = 0, any.missing = FALSE,         len = 1)    if (!is.null(scaleC)) {        checkmate::assertNumeric(scaleC, lower = 0, any.missing = FALSE)    }    checkmate::assertNumeric(scaleTo, len = 1, lower = 0, any.missing = FALSE)    checkmate::assertNumeric(gradTo, len = 1, lower = 0, any.missing = FALSE)    .ret <- list(covMethod = covMethod, typsize = typsize, fscale = fscale,         print.level = print.level, ndigit = ndigit, gradtol = gradtol,         stepmax = stepmax, steptol = steptol, iterlim = iterlim,         check.analyticals = check.analyticals, optExpression = optExpression,         literalFix = literalFix, literalFixRes = literalFixRes,         sumProd = sumProd, rxControl = rxControl, returnNlm = returnNlm,         stickyRecalcN = as.integer(stickyRecalcN), maxOdeRecalc = as.integer(maxOdeRecalc),         odeRecalcFactor = odeRecalcFactor, indTolRelax = indTolRelax,         eventType = eventType, shiErr = shiErr, shi21maxFD = as.integer(shi21maxFD),         optimHessType = optimHessType, hessErr = hessErr, shi21maxHess = as.integer(shi21maxHess),         censOption = censOption, eventSens = eventSens, sensMethod = sensMethod,         eventSens = eventSens, sensMethod = sensMethod, iterPrintControl = .iterPrintControl,         scaleType = scaleType, normType = normType, scaleCmax = scaleCmax,         scaleCmin = scaleCmin, scaleC = scaleC, scaleTo = scaleTo,         gradTo = gradTo, addProp = match.arg(addProp), calcTables = calcTables,         compress = compress, solveType = solveType, ci = ci,         sigdig = sigdig, sigdigTable = sigdigTable, genRxControl = .genRxControl,         boundedTransform = boundedTransform)    class(.ret) <- "nlmControl"    .ret})(covMethod = "nlm", typsize = NULL, fscale = 1, print.level = 0,     ndigit = 4, gradtol = 1e-06, stepmax = NULL, steptol = 1e-06,     iterlim = 10000, check.analyticals = FALSE, optExpression = TRUE,     literalFix = TRUE, literalFixRes = TRUE, sumProd = FALSE,     rxControl = structure(list(scale = NULL, method = c(liblsoda = 2L),         atol = 5e-07, rtol = 5e-07, maxsteps = 70000L, hmin = 0,         hmax = NA_real_, hini = 0, maxordn = 12L, maxords = 5L,         covsInterpolation = c(locf = 1L), addCov = TRUE, returnType = c(rxSolve = 0L),         sigma = NULL, sigmaDf = NULL, nCoresRV = 1L, sigmaIsChol = FALSE,         sigmaSeparation = "auto", sigmaXform = c(identity = 4L),         nDisplayProgress = 10000L, amountUnits = NA_character_,         timeUnits = "hours", addDosing = FALSE, stateTrim = Inf,         updateObject = FALSE, omega = NULL, omegaDf = NULL, omegaIsChol = FALSE,         omegaSeparation = "auto", omegaXform = c(variance = 6L),         nSub = 1L, thetaMat = NULL, thetaDf = NULL, thetaIsChol = FALSE,         nStud = 1L, dfSub = 0, dfObs = 0, seed = NULL, nsim = NULL,         minSS = 10L, maxSS = 10000L, strictSS = 1L, infSSstep = 12,         istateReset = TRUE, subsetNonmem = TRUE, hmaxSd = 0,         maxAtolRtolFactor = 0.1, from = NULL, to = NULL, by = NULL,         length.out = NULL, iCov = NULL, keep = NULL, keepF = character(0),         drop = NULL, warnDrop = TRUE, omegaLower = -Inf, omegaUpper = Inf,         sigmaLower = -Inf, sigmaUpper = Inf, thetaLower = -Inf,         thetaUpper = Inf, indLinPhiM = 0L, indLinPhiTol = 1e-07,         indLinMatExpType = c(expokit = 2L), indLinMatExpOrder = 6L,         idFactor = TRUE, mxhnil = 0L, hmxi = 0, warnIdSort = TRUE,         ssAtol = 5e-05, ssRtol = 5e-05, safeZero = 1L, sumType = c(pairwise = 1L),         prodType = c(`long double` = 1L), resample = NULL, resampleID = TRUE,         maxwhile = 100000L, cores = 0L, atolSens = 1.58113883008419e-06,         rtolSens = 1.58113883008419e-06, ssAtolSens = 0.000210848251714291,         ssRtolSens = 0.000210848251714291, simVariability = NA,         nLlikAlloc = NULL, useStdPow = 0L, naTimeHandle = c(ignore = 1L),         addlKeepsCov = FALSE, addlDropSs = TRUE, ssAtDoseTime = TRUE,         ss2cancelAllPending = FALSE, naInterpolation = c(locf = 1L),         keepInterpolation = c(na = 2L), safeLog = 1L, safePow = 1L,         ssSolved = TRUE, linCmtSensType = c(auto = 100L), linCmtSensH = 1e-04,         linCmtGillFtol = 0, linCmtGillK = 20L, linCmtGillStep = 4,         linCmtGillRtol = 1.49011611938477e-08, linCmtShiErr = 1.49011611938477e-08,         linCmtShiMax = 20L, linCmtScale = c(0, 0, 0, 0, 0, 0,         0), linCmtHcmt = 1L, linCmtHmeanI = c(geometric = 2L),         linCmtHmeanO = c(geometric = 2L), linCmtSuspect = 1e-06,         linCmtForwardMax = 2L, indOwnAlloc = -1L, maxExtra = 1000L,         tolFactor = NULL, serializeFile = NULL, dense = FALSE,         cvodeLinSolver = c(dense = 1L), stiff2 = 0L, autoSwitchMaxStiff = 10L,         autoSwitchMaxNonstiff = 3L, autoSwitchStiffFirst = 0L,         autoSwitchNonstifftol = 0.9, autoSwitchStifftol = 0.9,         autoSwitchDtfac = 2, autoSwitchSwitchMax = 5L, useLinCmt = TRUE,         file = NULL, chunkSize = NULL, parallel = 0L, .zeros = NULL), class = "rxControl"),     returnNlm = FALSE, stickyRecalcN = 4L, maxOdeRecalc = 5L,     odeRecalcFactor = 3.16227766016838, indTolRelax = TRUE, eventType = 2L,     shiErr = 6.05545445239334e-06, shi21maxFD = 20L, optimHessType = 2L,     hessErr = 6.05545445239334e-06, shi21maxHess = 20L, censOption = 0L,     eventSens = "jump", sensMethod = "default", eventSens = "jump",     sensMethod = "default", iterPrintControl = structure(list(        every = 1L, ncol = 4L, headerEvery = 10L, useColor = TRUE,         simple = FALSE), class = c("iterPrintControl", "list"    )), scaleType = 2L, normType = 1L, scaleCmax = 1e+05, scaleCmin = 1e-05,     scaleC = NULL, scaleTo = 1, gradTo = 1, addProp = "combined2",     calcTables = TRUE, compress = FALSE, solveType = 3L, ci = 0.95,     sigdig = 4, sigdigTable = 4, genRxControl = TRUE, boundedTransform = TRUE): formal argument "eventSens" matched by multiple actual arguments

print(fit2)
#> Error: object 'fit2' not found

# you can also get the nlm output with fit2$nlm

fit2$nlm
#> Error: object 'fit2' not found

# The nlm control has been modified slightly to include
# extra components and name the parameters
# }
```
