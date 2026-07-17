# foceiLik.R -- a general FOCE-family per-subject log-likelihood built from an
# rxode2 UI model (issue #414).  Exposes the FOCEi inner problem (the same engine
# advi/vae/fsaem/npag reuse) as a public load / run / unload lifecycle for
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
                             maxOdeRecalc = 5L, odeRecalcFactor = 10^0.5) {
  .interaction <- if (likelihood %in% c("foce", "focep")) 0L else 1L
  .foce <- if (identical(likelihood, "focep")) "foce+" else "nonmem"
  foceiControl(rxControl = rxControl, maxOuterIterations = 0L,
               maxInnerIterations = 0L, covMethod = "",
               interaction = .interaction, foce = .foce,
               sumProd = sumProd, optExpression = optExpression,
               literalFix = literalFix, addProp = addProp,
               calcTables = FALSE, compress = FALSE, eventSens = eventSens,
               indTolRelax = indTolRelax, maxOdeRecalc = maxOdeRecalc,
               odeRecalcFactor = odeRecalcFactor, print = 0L)
}

#' Load a general FOCE-family likelihood into memory
#'
#' Compiles the inner (FOCEi sensitivity) model from an rxode2 UI model,
#' preprocesses the data, and sets up the FOCEi inner problem in memory so that
#' individual log-likelihoods can be evaluated repeatedly (in parallel per
#' subject) at supplied etas without recompiling -- the setup used internally by
#' `est="advi"`, `est="vae"` and the f-SAEM fast kernel, exposed here for
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
#' @param ... Additional solving/model options passed to `.foceiLikControl`
#'   (e.g. `optExpression`, `addProp`, `eventSens`).
#' @return Invisibly, a handle list with the loaded system's dimensions:
#'   `initPar` (the estimation-scale parameter vector at the model's initial
#'   estimates, a ready `theta` for [foceiLikRun()]), `npars`, `ntheta`,
#'   `neta`, `nid`, `thetaNames`, `etaNames`, `idLvl` and `likelihood`.
#' @seealso [foceiLikRun()], [foceiLikUnload()]
#' @export
#' @author Matthew L. Fidler
foceiLikLoad <- function(object, data,
                         likelihood = c("focei", "focep", "foce"),
                         rxControl = rxode2::rxControl(), ...) {
  likelihood <- match.arg(likelihood)
  if (!is.null(nlmixr2global$foceiLikEnv)) {
    stop("a general likelihood system is already loaded; call foceiLikUnload() first",
         call. = FALSE)
  }
  .ui <- rxode2::rxUiDecompress(rxode2::assertRxUi(object))
  .control <- .foceiLikControl(likelihood, rxControl, ...)
  .control$est <- "focei"
  # Run the standard pre-process hooks (bounded transforms, covariates,
  # zero-omega, literal fix) so the inner problem matches a real focei fit's
  # parameterization; the hooks mutate .env0$ui/data/control in place.
  .env0 <- new.env(parent = emptyenv())
  .env0$ui <- .ui
  .env0$data <- data
  .env0$control <- .control
  .preProcessHooksRun(.env0, "focei")
  .ui <- rxode2::rxUiDecompress(.env0$ui)
  .data <- .env0$data
  .control <- .env0$control
  # advi/fsaem-style inner setup on the hooked ui
  .ui$control <- .control
  .env <- .ui$foceiOptEnv
  .env$ui <- .ui
  .env$est <- "focei"
  .env$table <- NULL
  .foceiPreProcessData(.data, .env, .ui, .control$rxControl)
  .env$control$est <- "focei"
  .env$control$printTop <- FALSE
  if (is.null(.env$control$nF)) .env$control$nF <- 0L
  .env$control$needOptimHess <- isTRUE(any(.ui$predDfFocei$distribution != "norm"))
  .env$aqn <- 0L; .env$qx <- double(0); .env$qw <- double(0); .env$qfirst <- FALSE
  .env$nAGQ <- 0L; .env$aqLow <- -Inf; .env$aqHi <- Inf; .env$nEstOmega <- 0L
  .neta <- length(.env$etaNames)
  .nid <- length(.env$idLvl)
  .env$etaMat <- matrix(0, .nid, .neta)
  .initPar <- as.numeric(foceiLikLoad_(.env))
  .iniDf <- .ui$iniDf
  .handle <- list(initPar = .initPar,
                  npars = length(.initPar),
                  ntheta = sum(!is.na(.iniDf$ntheta)),
                  neta = .neta,
                  nid = .nid,
                  thetaNames = .env$thetaNames,
                  etaNames = .env$etaNames,
                  idLvl = .env$idLvl,
                  likelihood = likelihood)
  nlmixr2global$foceiLikEnv <- .handle
  invisible(.handle)
}

#' Evaluate a loaded general FOCE-family likelihood at supplied etas
#'
#' Writes the population parameter vector into the loaded system and returns the
#' per-subject log-likelihood at the supplied etas, computed in parallel over
#' subjects.  Requires a system loaded by [foceiLikLoad()].
#'
#' @param theta The estimation-scale parameter vector (length `handle$npars`),
#'   matching the FOCEi optimizer parameterization: population thetas followed
#'   by the estimated Omega elements.  `handle$initPar` from [foceiLikLoad()] is
#'   a ready starting value.
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
