# foceiLik.R -- a general FOCE-family per-subject log-likelihood built from an
# rxode2 UI model (issue #414).  Exposes the FOCEi inner problem (the same engine
# vi/vae/npag reuse) as a public load / run / unload lifecycle for
# MCMC/SAMBA-style callers that need individual log-likelihoods, evaluated in
# parallel per id, outside a full fit.  Only one system can be loaded at a time.

#' A foceiControl carrying the requested inner likelihood + solving options.
#'
#' Mirrors `.adviInnerFoceiControl`: focei -> interaction=1; focep -> foce+
#' (interaction=0, residual variance at the live conditional eta); foce ->
#' nonmem (interaction=0, R frozen at eta=0).  `maxInnerIterations=0` means the
#' inner is evaluated at the supplied etas, never re-optimized.
#' @noRd
.foceiLikControl <- function(likelihood, rxControl,
                             sumProd = FALSE, optExpression = TRUE,
                             literalFix = FALSE, addProp = "combined2",
                             eventSens = "jump", indTolRelax = TRUE,
                             maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5,
                             scaleType = "nlmixr2", scaleTo = 1.0,
                             fallbackFD = FALSE) {
  .interaction <- if (likelihood %in% c("foce", "focep")) 0L else 1L
  .foce <- if (identical(likelihood, "focep")) "foce+" else "nonmem"
  foceiControl(rxControl = rxControl, maxOuterIterations = 0L,
               maxInnerIterations = 0L, covMethod = "",
               interaction = .interaction, foce = .foce,
               sumProd = sumProd, optExpression = optExpression,
               literalFix = literalFix, addProp = addProp,
               calcTables = FALSE, compress = FALSE, eventSens = eventSens,
               indTolRelax = indTolRelax, maxOdeRecalc = maxOdeRecalc,
               odeRecalcFactor = odeRecalcFactor, print = 0L,
               scaleType = scaleType, scaleTo = scaleTo,
               fallbackFD = fallbackFD)
}

#' Load a general FOCE-family likelihood into memory
#'
#' Compiles the inner (FOCEi sensitivity) model from an rxode2 UI model,
#' preprocesses the data, and sets up the FOCEi inner problem in memory so that
#' individual log-likelihoods can be evaluated repeatedly (in parallel per
#' subject) at supplied etas without recompiling -- the setup used internally by
#' `est="emvi"`/`est="fbvi"`, `est="vae"` and the f-SAEM fast kernel, exposed here for
#' MCMC/SAMBA-style callers (issue #414).
#'
#' Only one likelihood system may be loaded at a time; loading errors if one is
#' already loaded.  Use [foceiLikRun()] to evaluate and [foceiLikUnload()] to
#' free.
#'
#' @param object An `rxode2`/`nlmixr2` UI model (a model function or its
#'   compiled UI).
#' @param data The estimation data (a data frame with the usual nlmixr2
#'   columns).
#' @param likelihood The individual likelihood type: `"focei"` (FOCE with
#'   interaction), `"focep"` (FOCE+, interaction off with the residual variance
#'   at the conditional eta) or `"foce"` (NONMEM-style FOCE, residual variance
#'   frozen at eta=0).
#' @param rxControl An [rxode2::rxControl()] object for the ODE solving options.
#' @param scale The parameter scale [foceiLikRun()]'s `theta` is on:
#'   `"focei"` (default) is the FOCEi estimation scale (the historical
#'   behavior); `"natural"` pins the scaling to the identity
#'   (`scaleType="mult"`, `scaleTo=0`), so the population thetas in `theta`
#'   (and in the handle's `initPar`) are the natural-scale values directly
#'   comparable with `ui$iniDf$est` -- no unscaling step for external callers
#'   (#939).  With either choice the trailing omega-block entries of the
#'   parameter vector remain in the internal `diagXform` parameterization of
#'   `chol(Omega^-1)`; `"natural"` leaves them unscaled but does not change
#'   that parameterization.
#' @param combSens When `TRUE`, use the combined eta+theta sensitivity
#'   build (#958): the INNER model itself carries the theta-sensitivity
#'   columns (appended after the FOCEi block), so one ODE integration per
#'   subject serves the value, `d/d(eta)` AND `d/d(theta)` -- the fused
#'   batch entry `nlmixr2FoceiCondBatchThetaGrad` reads all three from a
#'   single solve (~2x cheaper per evaluation than the two-model path).
#'   Implies a theta-sensitivity request; the separate theta-sensitivity
#'   model is neither built nor compiled.
#' @param thetaSens When `TRUE`, also build and wire the theta-sensitivity
#'   model (`d(f)/d(theta)`, `d(V)/d(theta)` forward sensitivities for the
#'   estimated non-mu structural and residual-error thetas), the model the
#'   imp/advi engines use for their outer population gradient.  Off by
#'   default: building it costs compile and solve time a value-only caller
#'   should not pay (#939).  The handle's `thetaSensIdx` reports which
#'   `ntheta` indices carry sensitivities (`integer(0)` when none do, e.g.
#'   when every theta is mu-referenced).
#' @param est Estimation-method name the standard pre-process hooks run
#'   against (default `"focei"`).  The hooks read the named method's
#'   capability attributes (for example whether a bounded-transform
#'   reparameterization is wanted), so a caller preparing the problem for a
#'   different consumer can name it here and get that consumer's
#'   preprocessing guarantees instead of focei's (#939).  The inner engine
#'   itself is unaffected.
#' @param ... Additional solving/model options passed to `.foceiLikControl`
#'   (e.g. `optExpression`, `addProp`, `eventSens`).
#' @return Invisibly, a handle list with the loaded system's dimensions:
#'   `initPar` (the parameter vector at the model's initial estimates on the
#'   requested `scale`, a ready `theta` for [foceiLikRun()]), `npars`,
#'   `ntheta`, `neta`, `nid`, `thetaNames`, `etaNames`, `idLvl`,
#'   `likelihood`, `scale`, `thetaSens` (theta sensitivities wired, by
#'   either build), `combSens` (the #958 combined build) and
#'   `thetaSensIdx`.
#' @seealso [foceiLikRun()], [foceiLikUnload()]
#'
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- 1
#'     tv <- 3.45
#'     add.sd <- 0.7
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' # Set the likelihood up in memory once; only one may be loaded at a time
#' h <- foceiLikLoad(one.cmt, theo_sd, "focei")
#'
#' # The handle carries the dimensions and a ready starting parameter vector
#' h$nid
#' h$neta
#' h$initPar
#'
#' # Individual joint log-likelihood at eta = 0, one value per subject
#' eta <- matrix(0, h$nid, h$neta)
#' foceiLikRun(h$initPar, eta)
#'
#' # Free it when done (loading again before this errors)
#' foceiLikUnload()
#' }
#' @export
#' @author Matthew L. Fidler
foceiLikLoad <- function(object, data,
                         likelihood = c("focei", "focep", "foce"),
                         rxControl = rxode2::rxControl(),
                         scale = c("focei", "natural"),
                         thetaSens = FALSE,
                         combSens = FALSE,
                         est = "focei", ...) {
  likelihood <- match.arg(likelihood)
  scale <- match.arg(scale)
  checkmate::assertLogical(thetaSens, len = 1, any.missing = FALSE)
  checkmate::assertCharacter(est, len = 1, any.missing = FALSE, min.chars = 1)
  if (!is.null(nlmixr2global$foceiLikEnv)) {
    stop("a general likelihood system is already loaded; call foceiLikUnload() first",
         call. = FALSE)
  }
  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(object))
  if (identical(scale, "natural")) {
    # identity scale/unscale: scaleType="mult" with scaleTo=0 returns the
    # parameter unchanged in both directions (see unscalePar()/scalePar(),
    # src/inner.cpp), so the estimation scale IS the natural scale (#939)
    .control <- .foceiLikControl(likelihood, rxControl,
                                 scaleType = "mult", scaleTo = 0, ...)
  } else {
    .control <- .foceiLikControl(likelihood, rxControl, ...)
  }
  .control$est <- "focei"
  # Run the standard pre-process hooks (bounded transforms, covariates,
  # zero-omega, literal fix) so the inner problem matches the named consumer's
  # parameterization; the hooks mutate .env0$ui/data/control in place.  `est`
  # names the method whose capability attributes the hooks consult -- for
  # the default "focei" this matches a real focei fit.
  .env0 <- new.env(parent = emptyenv())
  .env0$ui <- .ui
  .env0$data <- data
  .env0$control <- .control
  .preProcessHooksRun(.env0, est)
  .ui <- rxode2::rxUiDecompress(.env0$ui)
  .data <- .env0$data
  .control <- .env0$control
  # theta-sensitivity request: computed on the HOOKED ui (bounded transforms
  # may have changed the parameterization).  thetaSensLoad makes
  # .foceiOptEnvLik build the model and foceiSetup_/vaeInnerSetup_ wire its
  # lhs offsets; impThetaSensIdx (0-based) names the estimated non-mu thetas
  # that get sensitivity columns -- exactly the .adviInnerSetup arrangement.
  .thetaSensIdx <- integer(0)
  if (isTRUE(thetaSens)) {
    .thetaSensIdx <- as.integer(.impmapEstTheta(.ui)$all)
    .control$thetaSensLoad <- TRUE
    .control$impThetaSensIdx <- .thetaSensIdx - 1L
  }
  # combined eta+theta sensitivity build (#958): the INNER model carries the
  # theta columns, so one solve serves value + d/d(eta) + d/d(theta); implies
  # a theta-sensitivity request
  if (isTRUE(combSens)) {
    if (!isTRUE(thetaSens)) {
      .thetaSensIdx <- as.integer(.impmapEstTheta(.ui)$all)
      .control$thetaSensLoad <- TRUE
      .control$impThetaSensIdx <- .thetaSensIdx - 1L
    }
    .control$combSens <- TRUE
  }
  # vi-style inner setup on the hooked ui
  .ui$control <- .control
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- "focei"
  .env$table <- NULL
  .foceiPreProcessData(.data, .env, .ui, .control$rxControl)
  .env$control$est <- "focei"
  if (isTRUE(thetaSens) || isTRUE(combSens)) {
    # foceiSetup_ reads thetaSensLoad/impThetaSensIdx from e$control (foceiO);
    # make sure both are present there (not only on the pre-build .control) so
    # op_focei wires the offsets -- same defensive re-set as .adviInnerSetup.
    .env$control$thetaSensLoad <- TRUE
    .env$control$impThetaSensIdx <- .thetaSensIdx - 1L
    if (isTRUE(combSens)) .env$control$combSens <- TRUE
  }
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  .env$control$needOptimHess <- isTRUE(any(.ui$predDfFocei$distribution != "norm"))
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .neta <- length(.env$etaNames)
  .nid <- length(.env$idLvl)
  .env$etaMat <- matrix(0, .nid, .neta)
  .initPar <- as.numeric(foceiLikLoad_(.env))
  # report what was actually wired, not what was asked for: the sensitivity
  # model build is a tryCatch(NULL) in .foceiOptEnvLik, so a request can fail
  # (and with no eligible thetas there is nothing to differentiate)
  # handle$thetaSens reports whether theta sensitivities are WIRED, however
  # they are carried: the separate theta-sensitivity model, or (#958) the
  # combined build whose inner model holds the columns
  .thetaSensBuilt <- (isTRUE(thetaSens) && !is.null(.env$model$thetaSens)) ||
    (isTRUE(combSens) && length(.thetaSensIdx) > 0L)
  if (isTRUE(thetaSens) && !.thetaSensBuilt && length(.thetaSensIdx) > 0L) {
    warning("the theta-sensitivity model could not be built; handle$thetaSens is FALSE",
            call. = FALSE)
  }
  .iniDf <- .ui$iniDf
  .handle <- list(initPar = .initPar,
                  npars = length(.initPar),
                  ntheta = sum(!is.na(.iniDf$ntheta)),
                  neta = .neta,
                  nid = .nid,
                  thetaNames = .env$thetaNames,
                  etaNames = .env$etaNames,
                  idLvl = .env$idLvl,
                  likelihood = likelihood,
                  scale = scale,
                  thetaSens = .thetaSensBuilt,
                  combSens = isTRUE(combSens),
                  thetaSensIdx = .thetaSensIdx)
  nlmixr2global$foceiLikEnv <- .handle
  invisible(.handle)
}

#' Evaluate a loaded general FOCE-family likelihood at supplied etas
#'
#' Writes the population parameter vector into the loaded system and returns the
#' per-subject log-likelihood at the supplied etas, computed in parallel over
#' subjects.  Requires a system loaded by [foceiLikLoad()].
#'
#' @param theta The parameter vector (length `handle$npars`) on the scale the
#'   system was loaded with: population thetas followed by the estimated Omega
#'   elements.  With `foceiLikLoad(scale="focei")` (the default) this is the
#'   FOCEi estimation scale; with `scale="natural"` the theta entries are the
#'   natural-scale values directly comparable with `ui$iniDf$est` (the omega
#'   entries stay in the internal `diagXform` parameterization either way).
#'   `handle$initPar` from [foceiLikLoad()] is a ready starting value.
#' @param eta A `nid` by `neta` matrix of random effects (one row per subject,
#'   in the loaded system's subject order).
#' @param type `"joint"` (default) returns the individual joint log density
#'   `log p(y_i, eta_i)`; `"cond"` returns the conditional data log-likelihood
#'   `log p(y_i | eta_i)` alone.  See Details.
#' @param cores Number of threads for the parallel per-subject evaluation.
#' @return A named numeric vector (length `nid`, named by subject id) of
#'   per-subject log-likelihoods.
#'
#' @details
#'
#' Both types are evaluated at the etas you supply, so both use each subject's
#' individual predictions; neither is a population (eta = 0) quantity.  They
#' differ only by the random-effect prior term:
#'
#' - `"cond"` is the conditional data log-likelihood `log p(y_i | eta_i)`, the
#'   observation contribution alone.
#' - `"joint"` is `log p(y_i, eta_i) = log p(y_i | eta_i) + log p(eta_i)`, which
#'   adds the Gaussian random-effect prior
#'   `log p(eta_i) = -0.5 eta_i' Omega^-1 eta_i + 0.5 log|Omega^-1| - neta/2 log(2 pi)`.
#'
#' So `"joint"` minus `"cond"` is exactly `log p(eta_i)`.  `"joint"` is the
#' default because it is the usual target for MCMC/SAMBA-style samplers: as a
#' function of `eta_i` it is the individual's posterior kernel, and it is the
#' quantity the FOCEi inner problem optimizes over the etas.  Use `"cond"` when
#' you supply the random-effect density yourself, or when you need the
#' observation contribution separately.
#'
#' The prior is built from the loaded system's own `Omega^-1` and its log
#' determinant -- the same Omega the inner likelihood uses -- so `"joint"` stays
#' internally consistent with the engine rather than with the nominal `ini()`
#' values (the two differ by a small amount through Omega's internal
#' `rxSymInv` representation).
#'
#' For Gaussian endpoints the observation contribution follows nlmixr2's
#' internal residual-likelihood convention, `-0.5 err^2/r - 0.5 log(r)`, which
#' omits the additive `-0.5 log(2 pi)` per observation; general log-likelihood
#' (`ll()`) endpoints contribute the user's log density as written.  The eta
#' prior above is fully normalized.  Both types are therefore proper log
#' densities up to a fixed per-observation constant that does not depend on
#' `theta` or `eta`, so likelihood ratios, and any sampler that uses them, are
#' unaffected.
#'
#' @seealso [foceiLikLoad()], [foceiLikUnload()]
#'
#' @examples
#'
#' \donttest{
#'
#' one.cmt <- function() {
#'   ini({
#'     tka <- 0.45
#'     tcl <- 1
#'     tv <- 3.45
#'     add.sd <- 0.7
#'     eta.ka ~ 0.6
#'     eta.cl ~ 0.3
#'     eta.v ~ 0.1
#'   })
#'   model({
#'     ka <- exp(tka + eta.ka)
#'     cl <- exp(tcl + eta.cl)
#'     v <- exp(tv + eta.v)
#'     linCmt() ~ add(add.sd)
#'   })
#' }
#'
#' h <- foceiLikLoad(one.cmt, theo_sd, "focei")
#'
#' eta <- matrix(0, h$nid, h$neta)
#'
#' # The individual joint log density log p(y_i, eta_i) (the default)
#' foceiLikRun(h$initPar, eta)
#'
#' # The conditional data log-likelihood log p(y_i | eta_i) alone; the two
#' # differ by the Gaussian eta prior
#' foceiLikRun(h$initPar, eta, type = "cond")
#'
#' # Non-zero etas
#' set.seed(42)
#' foceiLikRun(h$initPar, matrix(stats::rnorm(h$nid * h$neta, 0, 0.1), h$nid, h$neta))
#'
#' # A new population parameter vector needs no reload
#' theta <- h$initPar
#' theta[1] <- theta[1] + 0.1
#' foceiLikRun(theta, eta)
#'
#' foceiLikUnload()
#' }
#' @export
#' @author Matthew L. Fidler
foceiLikRun <- function(theta, eta, type = c("joint", "cond"),
                        cores = rxode2::getRxThreads()) {
  type <- match.arg(type)
  .h <- nlmixr2global$foceiLikEnv
  if (is.null(.h)) {
    stop("no general likelihood system loaded; call foceiLikLoad() first",
         call. = FALSE)
  }
  theta <- as.numeric(theta)
  if (length(theta) != .h$npars) {
    stop(sprintf("'theta' must have length %d (the loaded system's npars)", .h$npars),
         call. = FALSE)
  }
  eta <- as.matrix(eta)
  if (ncol(eta) != .h$neta) {
    stop(sprintf("'eta' must have %d columns (one per random effect)", .h$neta),
         call. = FALSE)
  }
  if (nrow(eta) != .h$nid) {
    stop(sprintf("'eta' must have %d rows (one per subject)", .h$nid),
         call. = FALSE)
  }
  .cores <- as.integer(cores)
  if (is.na(.cores) || .cores < 1L) .cores <- 1L
  foceiLikSetTheta_(theta)
  .retType <- if (identical(type, "cond")) 1L else 0L
  .ll <- foceiLikEval_(eta, .cores, .retType)
  stats::setNames(.ll, .h$idLvl)
}

#' Unload the general FOCE-family likelihood from memory
#'
#' Frees the FOCEi inner problem set up by [foceiLikLoad()].  A no-op (returns
#' `FALSE`) if nothing is loaded.
#'
#' @return Invisibly `TRUE` if a system was freed, `FALSE` if none was loaded.
#' @seealso [foceiLikLoad()], [foceiLikRun()]
#' @export
#' @author Matthew L. Fidler
foceiLikUnload <- function() {
  if (is.null(nlmixr2global$foceiLikEnv)) {
    return(invisible(FALSE))
  }
  foceiLikUnload_()
  nlmixr2global$foceiLikEnv <- NULL
  invisible(TRUE)
}

#' Iteration printing / parameter history over a loaded likelihood
#'
#' The standard nlmixr2est iteration table (the shared `scale.h` machinery
#' behind every estimation method's printout) exposed for an EXTERNAL
#' sampler driving [foceiLikLoad()]'s conditional likelihood -- e.g.
#' nlmixr2stan's `est="stan"`.  `foceiLikIterPrintStart()` prints the header
#' and arms the row entry (entry 8 of the FOCEi C API table,
#' `nlmixr2FoceiIterPrintRow`), which the sampler's gradient evaluations
#' call from C; every `every`-th call prints a row AND records it, so the
#' returned history holds exactly the printed rows.
#' `foceiLikIterPrintEnd()` prints the closing line and returns the history
#' as a `parHistData`-style data frame.
#'
#' The display vector is CALLER-defined, not the internal focei parameter
#' vector: pass natural-scale thetas plus (optionally) the sampler's current
#' actual omega entries, named `om.<eta>` for variances and
#' `cov.<eta1>.<eta2>` for covariances.  (The internal vector's omega tail
#' is `chol(Omega^-1)` in the `diagXform` parameterization, which a
#' Bayesian sampler neither uses nor understands, so it is NOT printed.)
#'
#' @param every print/record cadence in row-entry calls (0 arms nothing;
#'   the row entry then records nothing and returns immediately)
#' @param initPar display vector at the starting point (defines the width)
#' @param names column names, same length as `initPar`
#' @param iterPrintControl optional list applied via the shared
#'   `scaleApplyIterPrintControl` (e.g. `list(useColor=FALSE, printNcol=6L)`)
#' @param xform optional back-transform list (the `.iterPrintXParFromUi`
#'   shape) driving the `X` row
#' @return `foceiLikIterPrintStart()`: invisibly `NULL`;
#'   `foceiLikIterPrintEnd()`: the recorded history data frame (or `NULL`
#'   when printing was never started)
#' @export
#' @author Matthew L. Fidler
foceiLikIterPrintStart <- function(every, initPar, names,
                                   iterPrintControl = NULL, xform = NULL) {
  checkmate::assertIntegerish(every, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assertNumeric(initPar, min.len = 1, any.missing = FALSE)
  checkmate::assertCharacter(names, len = length(initPar))
  invisible(foceiLikIterPrintStart_(as.integer(every), as.double(initPar),
                                    names, iterPrintControl, xform))
}

#' @rdname foceiLikIterPrintStart
#' @export
foceiLikIterPrintEnd <- function() {
  foceiLikIterPrintEnd_()
}
