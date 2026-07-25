#' Control for advi (automatic differentiation variational inference) in nlmixr2
#'
#' Variational-inference NLME estimation in the style of Kucukelbir et al.
#' (2017): the latent variables are transformed to an unconstrained real
#' coordinate space, a Gaussian variational family is posited there, and the ELBO
#' is maximized by stochastic gradient ascent using the reparameterization trick.
#'
#' `est="advi"` is a VARIATION of that algorithm, not the published one, in two
#' respects.  First, by default (`pointEstimate=TRUE`) it is a variational-EM
#' HYBRID: the variational family covers only the per-subject etas, while the
#' population parameters (thetas, omega, residual error) are point estimates
#' updated by an M-step.  Published ADVI puts EVERY parameter -- including the
#' covariance parameters -- in the variational family, so a full-rank family
#' there yields a joint posterior covariance over the whole parameter vector.
#' `pointEstimate=FALSE` is the closer analogue, adding the population vector to
#' the variational posterior under flat priors, though the between-subject omega
#' is still parameterized as per-eta log-variances rather than freely.  Second,
#' the gradient of the log-joint comes from the FOCEi forward sensitivities
#' (inner per-subject eta gradient plus the outer population sensitivity
#' contraction) rather than from automatic differentiation, so "AD" in the name
#' describes the algorithm family, not the derivative source.  Results should
#' therefore not be read as reproducing a Stan `vb()` fit except on models where
#' the two objectives coincide.  The whole optimization loop runs in C++.
#'
#' @inheritParams saemControl
#' @inheritParams foceiControl
#'
#' @param seed Random seed for the ADVI optimization (reparameterization
#'   sampling); default 42.  The Monte-Carlo gradient is stochastic, so a fixed
#'   seed makes every fit reproducible.  Reparameterization noise is drawn from a
#'   counter-based stream keyed by the global iteration index, so a shorter run
#'   is a bit-for-bit prefix of a longer one and results are independent of the
#'   number of cores.
#' @param iters Total number of ADVI (stochastic gradient ascent) iterations.
#' @param nMc Number of Monte-Carlo samples used to approximate the ELBO
#'   gradient at each iteration (the paper's `M`; typically 1-10).
#' @param adviFamily Variational family in the unconstrained space.
#'   `"fullRank"` (default) uses a block full-rank Gaussian: a dense
#'   `neta x neta` Cholesky factor per subject plus a dense block over the
#'   population vector (mean-field across blocks).  `"meanField"` uses a fully
#'   factorized (diagonal) Gaussian.  Mean-field is faster but is known to
#'   underestimate marginal variances.
#' @param pointEstimate When `TRUE` (default) run a variational-EM hybrid: the
#'   variational posterior covers the per-subject etas only, and the population
#'   parameters (thetas / omega / residual error) are point estimates maximized
#'   by the ELBO gradient; output semantics match FOCEi/SAEM.  When `FALSE` run
#'   full Bayes: the variational posterior also covers the unconstrained
#'   population vector, with flat priors.
#'
#'   Two things about "flat" are worth being explicit about, because they define
#'   the prior rather than merely describe the implementation.  (1) A BOUNDED
#'   theta is fitted on its unconstrained scale, and the log-determinant of that
#'   constraining transform IS added to the full-Bayes objective, so the flat
#'   prior is flat on the NATURAL parameter, as in Stan.  It is deliberately NOT
#'   added when `pointEstimate=TRUE`: a maximum-likelihood estimate has to stay
#'   invariant to reparameterization, which is why Stan's own `optimize`
#'   defaults to `jacobian=0`.  (2) The between-subject variances are carried as
#'   per-eta LOG-variances and no Jacobian is applied to them, so the prior is
#'   flat on `log(omega)` -- the conventional weakly-informative choice for a
#'   scale parameter, but a choice, not an accident: it is not flat on `omega`.
#'
#'   These point estimates maximize the ELBO, NOT the likelihood, and the
#'   difference is not merely cosmetic.  Since
#'   `ELBO = log p(y|theta) - KL(q || p(eta|y,theta))`, any dependence of that KL
#'   on `theta` displaces the maximizer from the MLE -- "variational maximum
#'   likelihood".  For VARIANCE components the displacement has a known
#'   direction: a variational family that understates posterior spread makes the
#'   omega M-step, `Omega = mean_i(mu_i mu_i' + Sigma_i)`, inherit that
#'   understatement, so between-subject variability is biased DOWNWARD.  The bias
#'   is worst for `adviFamily="meanField"`, which cannot represent within-subject
#'   posterior correlation at all; `"fullRank"` can, which is why it is the
#'   default.  Structural (typical-value) parameters are far less affected.  If
#'   the between-subject variances are themselves the quantity of interest,
#'   prefer `"fullRank"` and cross-check against `est="focei"` or `est="saem"`.
#' @param optim Stochastic optimizer.  `"advi"` (default) uses the paper's
#'   adaptive step-size sequence (Eqs 10-11); `"adam"` uses Adam.
#' @param adaptEta When `TRUE` (default) adaptively choose the step-size scale
#'   `eta` by a short search over `etaCandidates` before the main loop; when
#'   `FALSE` use a fixed `eta` (the first `etaCandidates` entry).
#' @param etaCandidates Candidate step-size scales searched when `adaptEta` is
#'   `TRUE`.  The default is narrower and smaller-valued than the paper's
#'   `c(0.01, 0.1, 1, 10, 100)` because these gradients come from FOCEi
#'   sensitivities on the model's own scale rather than from AD through a Stan
#'   program, so the useful step sizes sit lower; the paper's grid can be passed
#'   verbatim if wanted.  Each candidate costs `min(iters, 75)` iterations
#'   (a diverging one aborts early and is cheap), so the search is a substantial
#'   share of a fit -- widen it deliberately.  When the search selects the
#'   largest or smallest candidate the grid may be the binding constraint, and
#'   the fit says so in `$runInfo`; `$etaScores` reports the per-candidate
#'   scores behind the choice.
#' @param perNoCor Fraction of the run over which a declared correlated `omega`
#'   block is held at zero correlation, letting the population variances settle
#'   before the correlations are estimated.  This is
#'   \code{\link{saemControl}()}'s `perNoCor` rule (0.75 there as well); it has no
#'   effect on a model with no declared off-diagonals.
#'
#'    The resolved absolute iteration is stored with the fit and
#'   reused by `resume=`, so a resumed run releases the correlations at the same
#'   global iteration a single long run would -- recomputing the fraction from the
#'   resumed call's `iters` would re-apply a hold the original run had passed.
#' @param tau Stabilizing constant `tau > 0` in the step-size denominator
#'   (paper Eq 10); the step-size is insensitive to it.
#' @param alpha Weighting `alpha` in (0, 1) of new vs old gradient information in
#'   the step-size memory recursion (paper Eq 11).
#' @param tol Convergence tolerance on the relative change in the ELBO: the loop
#'   stops early once the change falls below this.  `NULL` (default) derives it
#'   from `sigdig` as `10^(-sigdig)`, the same rule `saemControl()` and
#'   `foceiControl()` use for their optimizer tolerances, so it tightens with
#'   `sigdig` instead of staying pinned.  `0` disables early stopping (run all
#'   `iters`).  Because the per-iteration ELBO is an `nMc`-sample
#'   Monte-Carlo estimate and therefore noisy, the test compares the MEAN over
#'   the last `evalElbo` iterations against the mean over the window before it,
#'   rather than consecutive iterations.  The `adaptEta` step-size search never
#'   stops early on this criterion (its scorer reads a short run as divergence).
#' @param klWarmup Number of iterations of PRIOR TEMPERING (0, the default,
#'   disables it).  During the warm-up the population prior is inflated by a
#'   factor ramping geometrically from `temperInit` down to 1, which down-weights
#'   the prior term of the ELBO and keeps the per-subject variational posterior
#'   from collapsing before it is informative.  This is the variational analogue
#'   of \code{\link{saemControl}()}'s `perSa` simulated-annealing phase (true
#'   simulated annealing does not transfer: variational inference has no MCMC
#'   kernel to keep wide).
#'
#'   It CHANGES THE OBJECTIVE FUNCTION MID-RUN.  Early iterations maximize a
#'   tempered surrogate rather than the ELBO, so the convergence theory for the
#'   ELBO does not cover the warm-up, the reported ELBO trace is not comparable
#'   across the boundary, and the `tol` early-stopping test is suppressed until
#'   tempering ends.  The `adaptEta` step-size search also scores candidates on
#'   the untempered objective.  Off by default for those reasons; turn it on for
#'   a model where the variational scale collapses early.
#' @param temperInit Initial prior inflation factor for `klWarmup` tempering
#'   (default 10); ignored when `klWarmup = 0`.
#' @param evalElbo Window length, in iterations, for the `tol` convergence test.
#'   Stan's ADVI re-evaluates the ELBO every 100 iterations with fresh draws;
#'   averaging the draws already taken is the cheaper equivalent.  Shrunk
#'   automatically on a short run so `iters` well below `evalElbo` can still
#'   trigger the check.
#' @param likelihood Inner likelihood used for the per-subject objective and
#'   gradient, run through the FOCEi inner interface: `"focei"` (default),
#'   `"foce"`, `"focep"`, or `"laplace"`.
#' @param returnAdvi When `TRUE` return the raw ADVI optimization object instead
#'   of the nlmixr2 fit.
#' @param resume Optional warm-resume state: a previous `est="advi"` fit (or its
#'   `$env$adviState`).  The optimization continues from that state for `iters`
#'   more iterations, bit-for-bit identical to a single fresh run of the combined
#'   length (the counter-based RNG is keyed by the global iteration index).
#'
#'   That equivalence requires every schedule point to be an ABSOLUTE iteration.
#'   A FRACTIONAL `perNoCor` cannot provide it, and not because of any bookkeeping
#'   that could be fixed: `perNoCor = 0.75` of one 120-iteration run releases the
#'   correlations at iteration 90, while 0.75 of a first 60-iteration leg releases
#'   at 45.  Those are different schedules, and the resulting correlation estimates
#'   genuinely differ.  Pin the schedule (`perNoCor = 90`) whenever a fit may be
#'   resumed; the resolved value is then stored with the fit and reused.
#'
#' @return advi control structure (class `adviControl`)
#' @export
#' @author Matthew L. Fidler
adviControl <- function(seed = 42L,
                        iters = 300L,
                        nMc = 1L,
                        adviFamily = c("fullRank", "meanField"),
                        pointEstimate = TRUE,
                        optim = c("advi", "adam"),
                        adaptEta = TRUE,
                        perNoCor = 0.75,
                        etaCandidates = c(0.01, 0.025, 0.05, 0.1, 0.25),
                        tau = 1.0,
                        alpha = 0.1,
                        tol = NULL,
                        evalElbo = 100L,
                        klWarmup = 0L,
                        temperInit = 10,
                        likelihood = c("focei", "foce", "focep", "laplace"),
                        returnAdvi = FALSE,
                        resume = NULL,

                        print = 1L,
                        useColor = NULL,
                        printNcol = NULL,

                        covMethod = c("advi", "analytic", "r,s", "r", "s", ""),
                        optExpression = TRUE,
                        sumProd = FALSE,
                        literalFix = TRUE,
                        literalFixRes = TRUE,
                        addProp = c("combined2", "combined1"),
                        calcTables = TRUE,
                        compress = FALSE,
                        adjObf = TRUE,
                        ci = 0.95,
                        sigdig = 4,
                        sigdigTable = NULL,

                        stickyRecalcN = 4,
                        maxOdeRecalc = 5,
                        odeRecalcFactor = 10^(0.5),
                        indTolRelax = TRUE,
                        eventSens = c("jump", "fd"),
                        rxControl = NULL,
                        ...) {

  checkmate::assertIntegerish(seed, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(iters, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(nMc, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertLogical(pointEstimate, len = 1, any.missing = FALSE)
  checkmate::assertLogical(adaptEta, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(etaCandidates, lower = 0, finite = TRUE, any.missing = FALSE, min.len = 1)
  checkmate::assertNumeric(tau, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  ## no upper bound: <= 1 is a fraction of the phase, > 1 is an ABSOLUTE
  ## iteration count (the only form under which a resumed fit can reproduce a
  ## single run -- a fraction of a shorter leg is a different schedule)
  checkmate::assertNumeric(perNoCor, any.missing = FALSE, lower = 0, finite = TRUE, len = 1)
  checkmate::assertNumeric(alpha, lower = 0, upper = 1, any.missing = FALSE, len = 1)
  ## tol follows the package-wide sigdig convention (10^-sigdig), the same
  ## derivation saemControl()/foceiControl() use for their optimizer tolerances,
  ## so raising sigdig tightens the ELBO convergence test with everything else
  ## rather than leaving it pinned at a value that silently stops matching.
  ## sigdig's own assertion permits NA, which would make 10^-sigdig NA and fail
  ## the tol assertion below with an unhelpful message -- fall back instead.
  ## The derivation runs BEFORE sigdig's own assertion further down, so guard on
  ## the type here too: without it adviControl(sigdig = "bad") reports a base
  ## arithmetic error from 10^-sigdig rather than the intended checkmate message.
  if (is.null(tol)) {
    .sdOk <- !is.null(sigdig) && is.numeric(sigdig) && length(sigdig) == 1L &&
      !is.na(sigdig)
    tol <- if (.sdOk) .sigdigOptTol(sigdig) else 1e-4
  }
  checkmate::assertNumeric(tol, lower = 0, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(evalElbo, lower = 1, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(klWarmup, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(temperInit, lower = 1, finite = TRUE, any.missing = FALSE, len = 1)
  checkmate::assertLogical(returnAdvi, len = 1, any.missing = FALSE)
  checkmate::assertLogical(optExpression, len = 1, any.missing = FALSE)
  checkmate::assertLogical(sumProd, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFix, len = 1, any.missing = FALSE)
  checkmate::assertLogical(literalFixRes, len = 1, any.missing = FALSE)
  checkmate::assertLogical(calcTables, len = 1, any.missing = FALSE)
  checkmate::assertLogical(compress, len = 1, any.missing = TRUE)
  checkmate::assertLogical(adjObf, len = 1, any.missing = TRUE)
  checkmate::assertIntegerish(stickyRecalcN, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assertIntegerish(maxOdeRecalc, any.missing = FALSE, len = 1)
  checkmate::assertNumeric(odeRecalcFactor, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assertLogical(indTolRelax, len = 1, any.missing = FALSE)
  adviFamily <- match.arg(adviFamily)
  optim <- match.arg(optim)
  likelihood <- match.arg(likelihood)
  # match.arg cannot match ""; treat it (skip covariance) like foceiControl does
  if (length(covMethod) == 1L && covMethod == "") {
    covMethod <- ""
  } else {
    covMethod <- match.arg(covMethod)
  }
  addProp <- match.arg(addProp)
  eventSens <- match.arg(eventSens)

  .xtra <- list(...)
  .bad <- names(.xtra)
  .bad <- .bad[!(.bad %in% c("genRxControl", "iterPrintControl"))]
  if (length(.bad) > 0) {
    stop("unused argument: ", paste(paste0("'", .bad, "'"), collapse = ", "),
         call. = FALSE)
  }

  .genRxControl <- FALSE
  if (!is.null(.xtra$genRxControl)) {
    .genRxControl <- .xtra$genRxControl
  }
  if (is.null(rxControl)) {
    if (!is.null(sigdig)) {
      rxControl <- .rxControlScaleSigdig(rxode2::rxControl(sigdig = sigdig), sigdig)
    } else {
      rxControl <- rxode2::rxControl(atol = 1e-4, rtol = 1e-4)
    }
    .genRxControl <- TRUE
  } else if (inherits(rxControl, "rxControl")) {
  } else if (is.list(rxControl)) {
    rxControl <- .rxControlScaleSigdig(do.call(rxode2::rxControl, rxControl), sigdig, skip = names(rxControl))
  } else {
    stop("solving options 'rxControl' needs to be generated from 'rxode2::rxControl'",
         call. = FALSE)
  }
  if (!is.null(sigdig)) {
    checkmate::assertNumeric(sigdig, lower = 1, finite = TRUE, any.missing = TRUE, len = 1)
    if (is.null(sigdigTable)) {
      sigdigTable <- round(sigdig)
    }
  }
  if (is.null(sigdigTable)) {
    sigdigTable <- 3
  }
  checkmate::assertIntegerish(sigdigTable, lower = 1, len = 1, any.missing = FALSE)

  .iterPrintControl <- .absorbIterPrintControl(print = print,
                                               printNcol = printNcol,
                                               useColor = useColor,
                                               iterPrintControl = .xtra$iterPrintControl)

  .ret <- list(seed = as.integer(seed),
               iters = as.integer(iters),
               nMc = as.integer(nMc),
               adviFamily = adviFamily,
               pointEstimate = pointEstimate,
               optim = optim,
               adaptEta = adaptEta,
               perNoCor = perNoCor,
               etaCandidates = as.numeric(etaCandidates),
               tau = tau,
               alpha = alpha,
               tol = tol,
               evalElbo = as.integer(evalElbo),
               klWarmup = as.integer(klWarmup),
               temperInit = as.numeric(temperInit),
               likelihood = likelihood,
               returnAdvi = returnAdvi,
               resume = resume,
               covMethod = covMethod,
               optExpression = optExpression,
               sumProd = sumProd,
               literalFix = literalFix,
               literalFixRes = literalFixRes,
               addProp = addProp,
               calcTables = calcTables,
               compress = compress,
               adjObf = adjObf,
               ci = ci,
               sigdig = sigdig,
               sigdigTable = sigdigTable,
               stickyRecalcN = as.integer(stickyRecalcN),
               maxOdeRecalc = as.integer(maxOdeRecalc),
               odeRecalcFactor = odeRecalcFactor,
               indTolRelax = indTolRelax,
               eventSens = eventSens,
               iterPrintControl = .iterPrintControl,
               rxControl = rxControl,
               genRxControl = .genRxControl)
  class(.ret) <- "adviControl"
  .ret
}

#' @export
rxUiDeparse.adviControl <- function(object, var) {
  .default <- adviControl()
  object$resume <- NULL                     # not deparsable (may be a whole fit)
  .w <- .deparseDifferent(.default, object, "genRxControl")
  .deparseFinal(.default, object, .w, var)
}

#' @rdname nmObjHandleControlObject
#' @export
nmObjHandleControlObject.adviControl <- function(control, env) {
  assign("adviControl", control, envir = env)
}

#' @rdname nmObjGetControl
#' @export
nmObjGetControl.advi <- function(x, ...) {
  .env <- x[[1]]
  if (exists("adviControl", .env)) {
    .control <- get("adviControl", .env)
    if (inherits(.control, "adviControl")) return(.control)
  }
  if (exists("control", .env)) {
    .control <- get("control", .env)
    if (inherits(.control, "adviControl")) return(.control)
  }
  stop("cannot find advi related control object", call. = FALSE)
}

#' @rdname getValidNlmixrControl
#' @export
getValidNlmixrCtl.advi <- function(control) {
  .ctl <- control[[1]]
  if (is.null(.ctl)) .ctl <- adviControl()
  if (is.null(attr(.ctl, "class")) && is(.ctl, "list")) .ctl <- do.call("adviControl", .ctl)
  if (!inherits(.ctl, "adviControl")) {
    .minfo("invalid control for `est=\"advi\"`, using default")
    .ctl <- adviControl()
  } else {
    .ctl <- do.call(adviControl, .ctl)
  }
  .ctl
}

#' @rdname nlmixr2Est
#' @export
nlmixr2Est.advi <- function(env, ...) {
  .ui <- env$ui
  rxode2::assertRxUiRandomOnIdOnly(.ui, " for the estimation routine 'advi'", .var.name = .ui$modelName)
  ## advi has no mixture support (mean-field VI on the theta-sensitivity model);
  ## reject mix() up front instead of running a wrong fit that fails late in the
  ## output tables with a cryptic "probabilities in a mixture ... sum to 0".
  rxode2::assertRxUiNoMix(.ui, " for the estimation routine 'advi'", .var.name = .ui$modelName)
  ## absorb the validated control (set by getValidNlmixrControl before dispatch)
  if (exists("control", envir = env) && inherits(env$control, "adviControl")) {
    assign("adviControl", env$control, envir = env)
  } else {
    assign("adviControl", adviControl(), envir = env)
  }
  ## Seed the ENTIRE estimation ONCE here and restore the caller's global RNG
  ## state afterward; the counter-based reparameterization stream inside the loop
  ## makes a shorter run a bit-for-bit prefix of a longer one.
  rxode2::rxWithSeed(env$adviControl$seed, {
    .adviFitModel(env)
  })
}
attr(nlmixr2Est.advi, "covPresent") <- TRUE
## ADVI optimizes in the unconstrained real coordinate space
attr(nlmixr2Est.advi, "unbounded") <- TRUE
attr(nlmixr2Est.advi, "iov") <- TRUE
