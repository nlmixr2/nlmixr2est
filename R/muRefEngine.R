#' Restart-loop engine for the mu-referenced FOCEI family
#'
#' Implements the restart cycle described in the mu-referenced FOCEI family
#' design: (1) run the ordinary FOCEI inner+outer problem with mu-ref
#' covariate thetas held fixed (excluded from the outer optimizer), (2) fit
#' a quick linear model (`muModel="lin"`) or reweighted variant
#' (`muModel="irls"`) of the back-calculated per-subject `phi` on the
#' covariate(s) to get new population/covariate-coefficient thetas and a
#' residual eta, (3) plug the residual etas into a real (`maxOuterIterations
#' =0`, normal `maxInnerIterations`) optimizing pass to get the actual
#' conditional-mode phi/theta and objective function value, (4) warm-start
#' the next cycle from the *pre-regression* phi-derived (step 1) etas. Loops
#' until the objective function value converges (`muModelTol`) or
#' `muModelMaxCycles` is reached, then finalizes with the mu-ref covariate
#' thetas unlocked so the final Hessian/covariance/SE calculation is
#' standard FOCEI.
#'
#' Implementation note: each restart-cycle fit is a genuine, independent
#' top-level `nlmixr()` call (on a freshly-decompressed copy of the ui with
#' the mu-ref covariate thetas marked `fix=TRUE` via `iniDf`, exactly the
#' same mechanism `ini(... ~ fix)` uses), not a raw re-invocation of
#' `.foceiFitInternal()`/`foceiFitCpp_()` on a reused environment.
#' `.foceiFitInternal()`/`foceiFitCpp_()` mutate and populate output-only
#' fields directly onto the environment they are given (`etaObf`,
#' `objective`, internal C++ setup state, ...) and are only really designed
#' to be called once per environment (the existing "theta reset" retry loop
#' only re-invokes after a *failed* call, never after a successful one).
#' Going through the full `nlmixr()` pipeline for each cycle is slower but
#' reuses the same well-tested setup path every other estimation method
#' already relies on. Only the *final*, single fit (after the restart loop
#' converges) uses the original, never-yet-used `env` passed in here,
#' matching how every other method calls `.foceiFitInternal()` exactly
#' once.
#'
#' @param env the focei setup environment (`ui$foceiOptEnv`), same shape
#'   `.foceiFitInternal()` normally receives
#' @return the fitted result environment, same shape `.foceiFitInternal()`
#'   normally returns
#' @author Matthew L. Fidler
#' @noRd
.muRefFitInternal <- function(env) {
  .ui <- env[["ui"]]
  .groups <- .muRefGroups(.ui)
  if (length(.groups) == 0L) {
    # nothing mu-ref-covariate-eligible with an eta -- behave like ordinary
    # focei
    return(.foceiFitInternal(env))
  }
  .control <- env[["control"]]
  .thetaNames <- env[["thetaNames"]]
  .data <- env[["dataSav"]]

  .fixNames <- unique(c(
    vapply(.groups, function(g) g$theta, character(1)),
    unlist(lapply(.groups, function(g) g$covariates$covariateParameter))
  ))

  # user-fixed covariate coefficients (ini(... ~ fix)) are respected by the
  # linear-model step, not re-estimated -- keyed by covariate name (not
  # theta name), matching .muRefLin()/.muRefIrls()'s fixedCoef contract
  .iniDf0 <- .ui$iniDf
  .userFixedTheta <- setNames(.iniDf0$fix[!is.na(.iniDf0$ntheta)],
                               .iniDf0$name[!is.na(.iniDf0$ntheta)])

  .covByOwner <- .data[!duplicated(.data$ID), , drop=FALSE]
  .covByOwner <- .covByOwner[order(.covByOwner$ID), , drop=FALSE]

  .muModel <- .control$muModel
  .fitFn <- if (identical(.muModel, "irls")) .muRefIrls else .muRefLin

  .origThetaFixed <- env[["thetaFixed"]]

  # build a fresh, independently-decompressed copy of the ui for a cycle
  # fit: same model, iniDf$fix set for the mu-ref covariate thetas, iniDf
  # est set to the requested starting values. muModel="none" so this
  # doesn't recursively re-enter the restart loop.
  .cycleUi <- function(thetaIni) {
    .u <- rxode2::rxUiDecompress(rxode2::rxUiCompress(.ui))
    .id <- .u$iniDf
    .w <- !is.na(.id$ntheta)
    .id$est[.w] <- as.numeric(thetaIni[.id$name[.w]])
    .id$fix[.w] <- .id$name[.w] %in% .fixNames
    rxode2::ini(.u) <- .id
    .u
  }

  .runFit <- function(thetaIni, etaMat, maxOuterIterations) {
    .u <- .cycleUi(thetaIni)
    .res <- try(
      suppressWarnings(suppressMessages(
        nlmixr2est::nlmixr(.u, .data, est="focei",
                           control=foceiControl(etaMat=etaMat,
                                                 maxOuterIterations=maxOuterIterations,
                                                 maxInnerIterations=.control$maxInnerIterations,
                                                 muModel="none",
                                                 print=0))
      )),
      silent=TRUE
    )
    if (inherits(.res, "try-error")) {
      stop("mu-referenced FOCEI family restart cycle failed: ",
           attr(.res, "condition")$message, call.=FALSE)
    }
    .res
  }

  .phiForGroup <- function(g, theta, etaDf) {
    .covDf <- as.data.frame(.covByOwner[, g$covariates$covariate, drop=FALSE])
    names(.covDf) <- g$covariates$covariate
    .covEffect <- rep(0, nrow(.covDf))
    for (.k in seq_len(nrow(g$covariates))) {
      .covEffect <- .covEffect +
        .covDf[[g$covariates$covariate[.k]]] * theta[[g$covariates$covariateParameter[.k]]]
    }
    list(phi = theta[[g$theta]] + .covEffect + etaDf[[g$eta]], cov = .covDf)
  }

  .curThetaIni <- setNames(env[["thetaIni"]], .thetaNames)
  .etaMat <- env[["etaMat"]]
  .prevObj <- Inf
  .maxCycles <- .control$muModelMaxCycles
  .tol <- .control$muModelTol

  for (.cycle in seq_len(.maxCycles)) {
    # Step 1: standard inner+outer with mu-ref covariate thetas fixed
    .res1 <- .runFit(.curThetaIni, .etaMat, .control$maxOuterIterations)
    .theta1 <- .res1$theta
    .eta1 <- .res1$eta # data.frame: ID, <eta columns>
    .etaNames1 <- setdiff(names(.eta1), "ID")
    .step1EtaMat <- as.matrix(.eta1[, .etaNames1, drop = FALSE])
    rownames(.step1EtaMat) <- NULL

    # Step 2: quick linear model per group -> new theta_pop/theta_cov/eta
    .newThetaIni <- .curThetaIni
    .newEtaMat <- .step1EtaMat
    for (g in .groups) {
      .pc <- .phiForGroup(g, .theta1, .eta1)
      .fixedCoef <- NULL
      for (.k in seq_len(nrow(g$covariates))) {
        .cp <- g$covariates$covariateParameter[.k]
        if (isTRUE(.userFixedTheta[[.cp]])) {
          .fixedCoef[g$covariates$covariate[.k]] <- .theta1[[.cp]]
        }
      }
      .rres <- .fitFn(.pc$phi, .pc$cov, fixedCoef = .fixedCoef)
      .newThetaIni[g$theta] <- .rres$theta
      for (.k in seq_len(nrow(g$covariates))) {
        .cov <- g$covariates$covariate[.k]
        .cp <- g$covariates$covariateParameter[.k]
        .newThetaIni[.cp] <- unname(.rres$coef[.cov])
      }
      .newEtaMat[, g$eta] <- .rres$eta
    }

    # Step 3: maxOuterIterations=0, normal inner optimization on the
    # residual etas -- the actual optimized phi/theta and OFV
    .res3 <- .runFit(.newThetaIni, .newEtaMat, 0L)
    .curObj <- .res3$objective

    .converged <- is.finite(.prevObj) &&
      abs(.curObj - .prevObj) < .tol * max(1, abs(.prevObj))
    .prevObj <- .curObj

    # Step 4: warm-start next cycle from the pre-regression (step 1) etas,
    # with the newly regressed thetas
    .curThetaIni <- .newThetaIni
    .etaMat <- .step1EtaMat

    if (.converged) break
  }

  # Finalize: unlock the mu-ref covariate thetas and run one ordinary
  # (unmodified) fit -- on the pristine, never-yet-used env passed in --
  # so covariance/SE proceeds exactly as standard FOCEI. thetaIni must stay
  # unnamed/positional (matching env[["thetaNames"]] order) -- that's the
  # format .foceiFitInternal()/foceiSetupTheta_() expect.
  .muRefFinalize(env, unname(.curThetaIni[.thetaNames]), .origThetaFixed, .etaMat,
                 .control$maxOuterIterations)
}

#' Finalize a mu-referenced FOCEI family restart loop
#'
#' Unlocks the mu-ref covariate thetas (no longer excluded from the outer
#' optimizer) and runs one ordinary FOCEI fit so the final Hessian/
#' covariance/SE calculation is 100% standard FOCEI -- the resulting fit
#' object is structurally indistinguishable from a normal FOCEI fit.
#'
#' @param env the focei setup environment (never yet used for a fit call)
#' @param thetaIni converged theta values from the restart loop
#' @param thetaFixed the *original* (pre-restart-loop) fixed vector, i.e.
#'   only reflecting the user's own `ini(... ~ fix)` declarations
#' @param etaMat converged etaMat from the restart loop
#' @param maxOuterIterations the user's originally-requested
#'   `maxOuterIterations`
#' @return the fitted result environment
#' @author Matthew L. Fidler
#' @noRd
.muRefFinalize <- function(env, thetaIni, thetaFixed, etaMat, maxOuterIterations) {
  env[["thetaIni"]] <- thetaIni
  env[["thetaFixed"]] <- thetaFixed
  env[["etaMat"]] <- etaMat
  .ctl <- env[["control"]]
  .ctl$maxOuterIterations <- maxOuterIterations
  env[["control"]] <- .ctl
  .foceiFitInternal(env)
}
