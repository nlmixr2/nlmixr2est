# Analytic FOCEI observed-information covariance (covType="analytic"): the exact
# R-matrix from analytic 1st/2nd-order sensitivities (rxode2 .rxSens) with the
# 3rd-order tensor by Shi (2021) finite differences (keeps the augmented ODE at
# O(ndir^2), not O(ndir^3)).  Any solve/scope failure returns NULL and the caller
# falls back to the finite-difference Hessian.
#
# The sensitivity machinery runs over a "direction set": one direction per eta
# (ETA_i_) plus one direction per non-mu-referenced structural theta (THETA_j_,
# e.g. a covariate coefficient).  A mu-referenced theta reuses its eta's
# direction, so a fully mu-referenced model has ndir == neta.

#' Install the stashed analytic covariance as `fit$cov` (the native path fills
#' only the theta block).  No-op on the FD fallback (which foceiCalcR already
#' warned about).  `covFull=FALSE` (default) installs the structural-theta submatrix
#' (NONMEM-matched theta cov, backwards-compatible shape); `covFull=TRUE` installs
#' the full theta+sigma+Omega matrix (identical theta SEs) -- the assembly is always
#' full.
#' @param .ret focei fit environment
#' @noRd
.foceiInstallAnalyticCov <- function(.ret) {
  # only covMethod="r" installs the analytic R^-1; "r,s"/"s" keep the native
  # sandwich / S-matrix cov (which the analytic R already fed via covR).
  if (!identical(as.integer(rxode2::rxGetControl(.ret$ui, "covMethod", 2L)), 2L)) return(invisible())
  if (!exists(".analyticCov", envir = .ret, inherits = FALSE)) return(invisible())
  .cov <- get(".analyticCov", envir = .ret)
  if (!is.matrix(.cov) || !all(is.finite(.cov))) return(invisible())
  .full <- isTRUE(rxode2::rxGetControl(.ret$ui, "covFull", FALSE))
  if (!.full && exists(".analyticThetaNames", envir = .ret, inherits = FALSE)) {
    .th <- get(".analyticThetaNames", envir = .ret)          # structural cov-theta block only
    .th <- .th[.th %in% rownames(.cov)]
    if (length(.th) > 0L) .cov <- .cov[.th, .th, drop = FALSE]
  }
  # PD guard: an indefinite (near-boundary) inverse installs negative variances ->
  # NaN SEs.  Reject and keep the native/FD cov rather than a plausible-looking wrong one.
  .ev <- suppressWarnings(eigen(.cov, symmetric = TRUE, only.values = TRUE)$values)
  if (any(diag(.cov) <= 0) || !all(is.finite(.ev)) || min(.ev) <= 0) {
    warning("analytic covariance is not positive definite; keeping the finite-difference covariance",
            call. = FALSE)
    return(invisible())
  }
  .ret$cov <- .cov                       # analytic-tier cov already carries dimnames
  # covFull=TRUE swaps in a larger matrix than C++ foceiFinalizeTables saw, so its
  # condition numbers (computed from the theta-only native cov) are stale -- recompute.
  if (.full && exists("objDf", envir = .ret, inherits = FALSE)) {
    .od <- get("objDf", envir = .ret)
    if ("Condition#(Cov)" %in% names(.od)) .od[["Condition#(Cov)"]] <- max(.ev) / min(.ev)
    if ("Condition#(Cor)" %in% names(.od)) {
      .evc <- suppressWarnings(eigen(stats::cov2cor(.cov), symmetric = TRUE, only.values = TRUE)$values)
      .od[["Condition#(Cor)"]] <- max(.evc) / min(.evc)
    }
    assign("objDf", .od, envir = .ret)
  }
  invisible()
}

#' Scope-relevant structural thetas + sensitivity direction map: one direction per
#' eta (ETA_i_) plus one per non-mu-ref structural theta (THETA_j_); a mu-ref theta
#' reuses its eta's direction (so a fully mu-ref model has ndir == neta).  Sigma and
#' IOV SD thetas (`sgName`/`iovVars`) are not structural.  `NULL` if none remain.
#' @noRd
.foceiAnalyticDirections <- function(ini, thetaForEta, sgName, neta, iovVars = character(0)) {
  thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]
  thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
  thStructRows <- thRows[!.iniIsFixed(ini, thRows$name) & !(thRows$name %in% sgName) &
                           !(thRows$name %in% iovVars), , drop = FALSE]
  thStruct <- thStructRows$name
  nth <- length(thStruct)
  if (nth == 0L) return(NULL)
  etaDirs <- paste0("ETA_", seq_len(neta), "_")
  dirTh <- integer(nth); nonMuTheta <- character(0)
  for (p in seq_len(nth)) {
    .mt <- which(thetaForEta == thStruct[p])        # eta(s) this theta mu-references
    if (length(.mt) > 1L) return(NULL)              # shared by >1 eta -> summing routes NYI, fall back to FD
    k <- if (length(.mt) == 1L) .mt else NA_integer_
    if (!is.na(k)) {
      dirTh[p] <- k                                 # reuse eta k's direction (free)
    } else {
      nonMuTheta <- c(nonMuTheta, paste0("THETA_", thStructRows$ntheta[p], "_"))
      dirTh[p] <- neta + length(nonMuTheta)         # its own true-sensitivity direction
    }
  }
  dirs <- c(etaDirs, nonMuTheta)
  list(thStruct = thStruct, thStructRows = thStructRows, dirs = dirs,
       dirTh = dirTh, ndir = length(dirs), nth = nth)
}

#' Sum the per-subject analytic observed-information R (theta+sigma+Omega) over
#' subjects; `NULL` on any build/solve/non-finite failure.  `startedEnv` (when
#' given) is flagged `.analyticStarted` before the first solve so the C++ hook
#' skips its FD fallback rather than run it on the replaced global solve.
#' @noRd
.foceiAnalyticAssembleR <- function(ui, th, ebes, ids, data, Om, ef, neta, nth, nsg, omd,
                                    dirs, dirTh, ndir, startedEnv = NULL, solveTol = 1e-10,
                                    iovDirScale = NULL, etaScale = NULL, interaction = 1L) {
  am <- .foceiAnalyticAugModelDirs(ui, dirs)
  if (is.null(am) || am$ndir != ndir) return(NULL)
  np <- nth + nsg + omd$nom
  etav <- paste0("ETA_", seq_len(neta), "_")
  .foce <- identical(as.integer(interaction), 0L)     # FOCE re-solves EBEs to S_FOCE=0
  # IOV reparameterization to xi = w * eta (unit occasion eta -> variance w^2): the
  # augmented model is solved at the ACTUAL (Param A) EBEs, then the occasion-eta
  # sensitivities are rescaled by 1/w (iovDirScale) and their EBEs by w (etaScale)
  # so the Omega-variance branch sees a genuine variance-w^2 random effect.
  rescale <- !is.null(iovDirScale) && any(iovDirScale != 1)
  R <- matrix(0, np, np)
  Oi <- solve(Om)                                     # constant across subjects (invert once)
  # `data` (dataSav) keys on the etTrans integer code (1..N); `ids` (etaObf$ID /
  # fit$eta$ID) is a factor whose LABELS are the original dataset IDs.  Key BOTH on
  # the integer code so non-1..N IDs (e.g. 101,102,..) and permutations join right.
  .byId <- split(data, as.character(data$ID))         # pre-split once (avoid per-subject rescan)
  .idCode <- if (is.factor(ids)) as.integer(ids) else match(ids, sort(unique(ids)))
  if (!is.null(startedEnv)) assign(".analyticStarted", TRUE, startedEnv)
  for (i in seq_along(ids)) {
    s <- .byId[[as.character(.idCode[i])]]             # subject's actual events (integer-code join)
    if (is.null(s) || nrow(s) == 0L) return(NULL)      # unmatched subject -> caller falls back to FD
    obs <- s[s$EVID == 0, , drop = FALSE]
    eta0 <- ebes[i, ]
    # FOCE (interaction=0): nlmixr's stored EBEs do NOT stationarize S_FOCE (the fit
    # inner problem is a separate nlmixr issue), so re-solve each subject to the FOCE
    # stationary point before forming R.  No-op for additive/FOCEI (|S_FOCE|~0 -> skip).
    if (.foce) {
      eta0 <- .foceiAnalyticFoceEbe(am, th, eta0, s, obs$TIME, obs$DV, etav, ef, Oi, neta, solveTol)
      if (is.null(eta0)) return(NULL)
    }
    p <- c(th, setNames(eta0, etav))                  # solve at the Param A (unit-occ-eta) EBEs
    E <- .foceiAnalyticSolveSubjectFD3(am, p, s, obs$TIME, tol = solveTol)
    if (is.null(E)) return(NULL)                       # solve failure -> caller falls back to FD
    E$y <- obs$DV
    ehat <- eta0
    if (rescale) {
      E$a <- sweep(E$a, 2, iovDirScale, `*`)           # a_B = a_A / w on occasion directions
      for (d1 in seq_len(ndir)) for (d2 in seq_len(ndir))
        E$A[, d1, d2] <- E$A[, d1, d2] * iovDirScale[d1] * iovDirScale[d2]
      for (d1 in seq_len(ndir)) for (d2 in seq_len(ndir)) for (d3 in seq_len(ndir))
        E$Ath[, d1, d2, d3] <- E$Ath[, d1, d2, d3] * iovDirScale[d1] * iovDirScale[d2] * iovDirScale[d3]
      ehat <- ehat * etaScale                          # xi_hat = w * eta_hat on occasion etas
    }
    Ri <- tryCatch(.foceiAnalyticSubjectR(E, ehat, Om, ef, neta, nth, nsg, ef$sgVar, omd,
                                          ndir = ndir, dirTh = dirTh, Oi = Oi, interaction = interaction),
                   error = function(e) NULL)
    if (is.null(Ri) || !all(is.finite(Ri))) return(NULL)
    R <- R + Ri
  }
  R
}

#' C++ `foceiCalcR` hook for covType="analytic": from the live focei env `e`,
#' assemble the full theta+sigma+Omega observed-information R, stash the natural
#' cov (installed by [.foceiInstallAnalyticCov]), and return the theta `R.0` for
#' the native SE path -- or `NULL` to fall back to the FD Hessian.
#' @noRd
.foceiCalcRanalytic <- function(e) {
  tryCatch({
    ui <- get("ui", e)
    # covType="analytic" opt-in; anything else keeps the finite-difference Hessian
    if (!identical(rxode2::rxGetControl(ui, "covType", "fd"), "analytic")) return(NULL)
    if (!.hasRxSens()) return(NULL)
    # bounded thetas are estimated on a transformed scale and the pre-final Jacobian
    # hook corrects env$cov to the natural scale; an analytic install would overwrite
    # that with the internal-scale cov -> bow out to the (Jacobian-correct) FD path.
    if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L) return(NULL)
    # scope: the analytic observed information is FOCEI-only (conditional, with the
    # interaction term) and a single Gaussian endpoint; anything else -> FD.
    if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(NULL)             # FO/FOI
    interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))                   # 1 FOCEI / 0 FOCE
    if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L) return(NULL)                  # AGQ quadrature
    ef <- .foceiAnalyticErrFull(ui)
    if (is.null(ef)) return(NULL)

    ini <- ui$iniDf
    .map <- .foceiEtaThetaMap(ui)
    etaNames <- .map$etaNames; neta <- length(etaNames)
    if (neta == 0L) return(NULL)
    thetaForEta <- .map$thetaForEta

    # IOV (inter-occasion variability): nlmixr rewrites `v ~ .. | occ` BEFORE the fit
    # into a SD-scale THETA `v` (in skipCov) plus k FIXED unit-variance per-occasion
    # etas `rx.v.<occ>` entering the predictor as `abs(v) * sum(rx.v.<occ>*(occ==j))`.
    # We treat the shared variance omega_v = v^2 as ONE Omega-variance parameter over
    # those k occasion etas (reusing the "om" branch), so the occasion etas are the
    # only permitted non-mu-referenced etas.
    iovVars <- .uiIovEnv$iovVars
    # only the SD-scale IOV predictor (abs(v)*occ-eta) is implemented; var/logsd/logvar
    # use a different predictor + chain rule -> bow out to FD.
    if (length(iovVars) > 0L &&
          !identical(rxode2::rxGetControl(ui, "iovXform", "sd"), "sd")) return(NULL)
    occEta <- lapply(iovVars, function(v) grep(paste0("^rx\\.", gsub(".", "\\.", v, fixed = TRUE), "\\."), etaNames))
    occEtaAll <- sort(unique(unlist(occEta)))
    naEta <- which(is.na(thetaForEta))
    # a plain non-mu-ref eta keeps its own ETA_i_ direction + eta-named Omega variance;
    # occasion (IOV) etas are the grouped shared-variance case handled below.
    if (length(occEtaAll) > 0L && !all(occEtaAll %in% naEta)) return(NULL)  # sanity: occ etas must be non-mu-ref

    # converged estimates (e$theta, NOT ui$iniDf$est which holds the initials)
    thNames <- get("thetaNames", e)
    thVals  <- get("theta", e)$theta
    names(thVals) <- thNames
    # cov parameters = non-skipped thetas (skipCov drops residual/fixed/IOV)
    skip <- get("skipCov", e)[seq_along(thNames)]
    covParams <- thNames[!skip]

    # converged residual sigma into the error-model evaluator (ef reads initials)
    .valc <- setNames(as.numeric(thVals[ef$sgName]), ef$sgVar)
    ef$ev <- local({ v <- .valc; function(expr, f, y) eval(expr, c(list(f = f, y = y), as.list(v))) })

    # scope: fixed structural thetas break the eta<->theta indexing
    if (any(.iniIsFixed(ini, thetaForEta))) return(NULL)
    keep <- !.iniIsFixed(ini, ef$sgName); ef$sgVar <- ef$sgVar[keep]; ef$sgName <- ef$sgName[keep]  # drop fixed sigma

    Om <- get("omega", e)
    # Param B: give each occasion-eta group its shared variance w^2 (the internal
    # model fixes them to unit variance and carries w=|v| in the predictor); the
    # occasion-eta sensitivities/EBEs are rescaled to match in .foceiAnalyticAssembleR.
    iovW <- vapply(iovVars, function(v) abs(as.numeric(thVals[v])), numeric(1))
    etaScale <- rep(1, neta)
    iovGroups <- vector("list", length(iovVars))
    for (g in seq_along(iovVars)) {
      w <- iovW[g]; oe <- occEta[[g]]
      if (length(oe) == 0L || !is.finite(w) || w <= 0) return(NULL)   # unidentified/zero IOV -> FD
      for (j in oe) Om[j, j] <- w^2
      etaScale[oe] <- w
      iovGroups[[g]] <- list(idx = oe, w = w)
    }
    pairs <- .foceiOmegaPairs(Om, ini)                # free Omega lower-triangle (declared blocks)
    # occasion-eta diagonals are grouped into the IOV param(s), so drop them from the
    # ordinary element list (they are FIXED and thus already dropped by .omegaFixed).
    if (length(occEtaAll) > 0L && nrow(pairs) > 0L)
      pairs <- pairs[!(pairs[, 1] %in% occEtaAll | pairs[, 2] %in% occEtaAll), , drop = FALSE]
    omd <- .omegaVarCovDeriv(Om, pairs, iov = if (length(iovGroups)) iovGroups else NULL)

    # Uniform direction assembly (shared with the standalone oracle).  IOV variance
    # thetas are NOT structural (predictor) thetas: they never get a THETA_j_ direction
    # (the occasion-eta sensitivities already carry the w factor).
    .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta, iovVars)
    if (is.null(.dir)) return(NULL)
    thStruct <- .dir$thStruct
    dirs <- .dir$dirs; dirTh <- .dir$dirTh; ndir <- .dir$ndir; nth <- .dir$nth
    # occasion-eta directions carry a[,occ] = w * a_B; rescale by 1/w to the variance-w^2 basis
    iovDirScale <- rep(1, ndir)
    for (g in seq_along(iovVars)) iovDirScale[occEta[[g]]] <- 1 / iovW[g]

    nsg <- length(ef$sgVar)
    # om-param order matches .omegaVarCovDeriv: ordinary Omega elements THEN the IOV
    # shared-variance params (reported on the SD scale as `v`, matching the theta name).
    onm <- ifelse(is.na(thetaForEta), etaNames, thetaForEta)   # orphan (non-mu-ref) eta named by the eta
    fullNm <- c(thStruct, ef$sgName, .foceiOmegaCovNames(pairs, onm), iovVars)  # full natural-scale order

    # cov-param order (op_focei) -> theta-block position (thStruct order)
    idx <- match(covParams, thStruct)
    if (anyNA(idx)) return(NULL)

    th <- setNames(as.numeric(thVals[thNames]), paste0("THETA_", seq_along(thNames), "_"))
    etav <- paste0("ETA_", seq_len(neta), "_")
    etaObf <- get("etaObf", e)
    ebes <- as.matrix(etaObf[, paste0("ETA[", seq_len(neta), "]"), drop = FALSE])
    ids  <- etaObf$ID
    data <- get("dataSav", e)
    if (!is.null(data$CENS) && any(data$CENS != 0, na.rm = TRUE)) return(NULL)  # censored -> FD

    # Full natural-scale observed-information R (theta + sigma + Omega), summed over
    # subjects.  startedEnv=e flags `.analyticStarted` before the augmented solve so
    # the C++ hook skips its finite-difference fallback if the assembly then fails.
    Rfull <- .foceiAnalyticAssembleR(ui, th, ebes, ids, data, Om, ef, neta, nth, nsg, omd,
                                     dirs = dirs, dirTh = dirTh, ndir = ndir,
                                     startedEnv = e, solveTol = .foceiAnalyticSolveTol(ui),
                                     iovDirScale = iovDirScale, etaScale = etaScale,
                                     interaction = interaction)
    if (is.null(Rfull)) return(NULL)
    dimnames(Rfull) <- list(fullNm, fullNm)

    # full natural-scale cov -> fit$cov; theta SEs flow through the native path
    .covNat <- tryCatch(solve(Rfull), error = function(e) NULL)
    if (is.null(.covNat)) return(NULL)
    dimnames(.covNat) <- list(fullNm, fullNm)
    # the theta-block R (op_focei cov-param order) decides success; the native
    # covR = Rinv path (issue #666 fix) reports these SEs directly.
    .ct <- .covNat[seq_len(nth), seq_len(nth), drop = FALSE][idx, idx, drop = FALSE]
    .R0 <- tryCatch(solve(.ct), error = function(e) NULL)
    if (is.null(.R0)) return(NULL)   # inversion failed -> FD, and do NOT leave a stale .analyticCov
    assign(".analyticCov", .covNat, envir = e)           # stash only after the deciding inversion
    # covFull=FALSE installs only the native cov-theta block: respect skipCov so the
    # shape matches covType="fd" (do NOT widen to skipped structural thetas).
    assign(".analyticThetaNames", thStruct[thStruct %in% covParams], envir = e)
    .R0
  }, error = function(e) NULL)
}

#' Is the rxode2 symbolic-sensitivity machinery available (scope gate)?
#' @return logical, TRUE when `rxExpandSens2_` and symengine are present
#' @noRd
.hasRxSens <- function() {
  .ns <- asNamespace("rxode2")
  exists("rxExpandSens2_", envir = .ns, inherits = FALSE) &&
    exists("rxOmegaVarCovDeriv", envir = .ns, inherits = FALSE) &&
    requireNamespace("symengine", quietly = TRUE)
}

#' Augmented-solve tolerance for the SEs: `covSolveTol` if the user set it, else
#' tightened from `sigdig` (default 1e-10).
#' @noRd
.foceiAnalyticSolveTol <- function(ui) {
  .user <- tryCatch(rxode2::rxGetControl(ui, "covSolveTol", NULL), error = function(e) NULL)
  if (!is.null(.user) && is.finite(.user) && .user > 0) return(.user)
  .sd <- suppressWarnings(as.numeric(rxode2::rxGetControl(ui, "sigdig", 4)))
  if (!is.finite(.sd)) .sd <- 4
  max(1e-14, min(1e-8, 10^-(.sd + 6)))
}

#' Non-Cholesky Omega derivatives for the analytic Omega block: Omega^{-1} and
#' log|Omega| first/second derivatives w.r.t. the free variance-covariance
#' elements `pairs` (each row `c(a, b)`, `a >= b`), from rxode2's
#' `rxOmegaVarCovDeriv`.  Returns `nom`, `dOi` (list of dOmega^{-1}/dw), `d2Oi`
#' (list-of-lists), and `d2LD` (the log-determinant second derivatives).
#'
#' `iov` (optional) is a list of IOV shared-variance parameters; each element is
#' `list(idx = <occasion-eta Omega positions>, w = <SD-scale value>)`.  The k
#' occasion diagonals all equal ONE parameter (variance omega = w^2), so its
#' Omega derivative is the SUM of the per-element derivatives over those positions
#' (the block is diagonal, so within-group cross second derivatives are summed and
#' cross-terms to the disjoint ordinary blocks vanish).  The group is chain-ruled
#' from the variance omega=w^2 to the reported SD scale w (`dOmega/dw = 2w`,
#' `d2Omega/dw^2 = 2`), so the resulting observed-information entry is directly
#' `var(w)` -- no separate delta method.  IOV params are appended AFTER `pairs`.
#' @noRd
.omegaVarCovDeriv <- function(Om, pairs, iov = NULL) {
  d <- rxode2::rxOmegaVarCovDeriv(Om, order = 2L)
  key <- function(P) paste(P[, 1], P[, 2], sep = "-")
  ek <- key(d$elements)
  ordIdx <- match(key(pairs), ek)
  nOrd <- length(ordIdx)
  nIov <- length(iov)
  nom <- nOrd + nIov
  Z <- matrix(0, nrow(Om), ncol(Om))
  dOi <- vector("list", nom)
  d2Oi <- lapply(seq_len(nom), function(i) lapply(seq_len(nom), function(j) Z))
  d2LD <- matrix(0, nom, nom)
  # ordinary variance/covariance parameters (one per row of `pairs`)
  if (nOrd > 0L) {
    for (a in seq_len(nOrd)) {
      dOi[[a]] <- d$dOmegaInv[[ordIdx[a]]]
      for (b in seq_len(nOrd)) d2Oi[[a]][[b]] <- d$d2OmegaInv[[ordIdx[a]]][[ordIdx[b]]]
    }
    d2LD[seq_len(nOrd), seq_len(nOrd)] <- d$d2LogDet[ordIdx, ordIdx, drop = FALSE]
  }
  # IOV shared-variance parameters (grouped + chain-ruled to the SD scale w)
  if (nIov > 0L) {
    for (g in seq_len(nIov)) {
      p <- nOrd + g
      w <- iov[[g]]$w
      gi <- match(key(cbind(iov[[g]]$idx, iov[[g]]$idx)), ek)   # (occ, occ) element rows
      sdOi <- Reduce(`+`, d$dOmegaInv[gi])                      # sum_e dOmegaInv_e
      dOi[[p]] <- 2 * w * sdOi                                  # dOmegaInv/dw
      s2 <- Z
      for (a in gi) for (b in gi) s2 <- s2 + d$d2OmegaInv[[a]][[b]]
      d2Oi[[p]][[p]] <- 4 * w^2 * s2 + 2 * sdOi                 # d2OmegaInv/dw^2
      d2LD[p, p] <- 4 * w^2 * sum(d$d2LogDet[gi, gi]) + 2 * sum(d$dLogDet[gi])
    }
  }
  list(nom = nom, dOi = dOi, d2Oi = d2Oi, d2LD = d2LD)
}

#' Symbolic error machinery (rho/p f-derivatives + sigma partials) for one
#' Gaussian add/prop endpoint; `NULL` for anything else (multi-endpoint, lnorm,
#' propT/propF, a DV transform, combined1, or pure proportional) -> FD fallback.
#' @noRd
.foceiAnalyticErrFull <- function(ui) {
  # single Gaussian endpoint only: multiple endpoints pool error rows against one
  # rx_pred_ (the wrong likelihood)
  if (!is.null(ui$predDf) && nrow(ui$predDf) != 1L) return(NULL)
  ini <- ui$iniDf; er <- ini[!is.na(ini$err), , drop = FALSE]
  # whitelist: only additive / proportional error; any other rider (lnorm,
  # propT/propF, or a DV/variance transform boxCox / yeoJohnson / pow) is out of scope
  if (!all(er$err %in% c("add", "prop"))) return(NULL)
  # model-declared addProp wins; the control applies only when the model says "default"
  addPr <- as.character(ui$predDf$addProp)
  if (length(addPr) != 1L || is.na(addPr) || addPr == "default") {
    addPr <- tryCatch(rxode2::rxGetControl(ui, "addProp", "combined2"), error = function(e) "combined2")
  }
  if (identical(addPr, "combined1")) return(NULL)
  addN <- er$name[er$err == "add"]; propN <- er$name[er$err == "prop"]
  hasA <- length(addN) == 1L; hasP <- length(propN) == 1L
  if (!hasA && !hasP) return(NULL)
  # pure proportional error has variance sp^2 f^2 -> 0 as f -> 0, so the observed
  # information (1/R terms) blows up near zero predictions; leave it to the FD cov.
  if (hasP && !hasA) return(NULL)
  Rstr <- if (hasA && hasP) "sa^2+sp^2*f^2" else if (hasP) "sp^2*f^2" else "sa^2"
  Rq <- parse(text = Rstr)[[1]]
  rhoE <- bquote(0.5 * ((y - f)^2 / .(Rq) + log(.(Rq))))
  pE <- bquote(1 / .(Rq) + 0.5 * (.(D(Rq, "f")) / .(Rq))^2)
  # FOCE (interaction=0) pieces: the inner problem drops dR/deta, so its gradient
  # coefficient is the least-squares part only, qE = drho/df with R frozen = -(y-f)/R
  # (q'=dq/df, q''=d2q/df2 build the FOCE inner Hessian Hf and its 3-tensor); the FOCE
  # Laplace determinant curvature is pFE = 1/R (no 0.5*(R'/R)^2 interaction term).
  qE <- bquote(-(y - f) / .(Rq))
  pFE <- bquote(1 / .(Rq))
  DD <- function(e, ...) { for (v in c(...)) e <- D(e, v); e }
  sgVar <- c(if (hasA) "sa", if (hasP) "sp"); sgName <- c(if (hasA) addN, if (hasP) propN)
  val <- setNames(c(if (hasA) er$est[er$name == addN], if (hasP) er$est[er$name == propN]), sgVar)
  sc <- list(r1 = DD(rhoE, "f"), r2 = DD(rhoE, "f", "f"), r3 = DD(rhoE, "f", "f", "f"),
             p = pE, p1 = DD(pE, "f"), p2 = DD(pE, "f", "f"),
             q0 = qE, q1 = DD(qE, "f"), q2 = DD(qE, "f", "f"),
             pF = pFE, pF1 = DD(pFE, "f"), pF2 = DD(pFE, "f", "f"))
  per <- list(); for (s in sgVar) per[[s]] <- list(rf = DD(rhoE, "f", s), rff = DD(rhoE, "f", "f", s),
    ps = DD(pE, s), pf = DD(pE, "f", s),
    qs = DD(qE, s), qsf = DD(qE, "f", s), psF = DD(pFE, s), pfF = DD(pFE, "f", s))
  pair <- list(); for (i in seq_along(sgVar)) for (j in i:length(sgVar)) { a <- sgVar[i]; b <- sgVar[j]
    pair[[paste0(a, b)]] <- list(rss = DD(rhoE, a, b), rfss = DD(rhoE, "f", a, b), pss = DD(pE, a, b),
      qss = DD(qE, a, b), pssF = DD(pFE, a, b)) }
  list(sgVar = sgVar, sgName = sgName, sc = sc, per = per, pair = pair,
       ev = function(e, f, y) eval(e, c(list(f = f, y = y), as.list(val))))
}

#' Augmented rxode2 model with state + 1st/2nd-order sensitivities and the
#' prediction chain f1/f2 over an arbitrary direction set `dirs` (each ETA_i_ or
#' THETA_j_).  State sensitivities use rxode2's `.rxSens`; the higher-order
#' prediction chain f1/f2 is built here.  O(ndir^2), so it compiles far more
#' cheaply than a 3rd-order model.
#' @return list(augMod, dirs, ndir, st, P2) or `NULL` on failure
#' @noRd
.foceiAnalyticAugModelDirs <- function(ui, dirs) {
  tryCatch({
    .s <- ui$loadPruneSens
    .st <- rxode2::rxStateOde(.s)
    rxode2::.rxJacobian(.s, c(.st, dirs))
    # 1st-order sensitivities are already expanded for the gradient (free); if
    # unavailable the model is not differentiable, so bail before the 2nd-order build.
    .s1 <- rxode2::.rxSens(.s, dirs)
    if (length(.s1) == 0L) return(NULL)
    .s2 <- rxode2::.rxSens(.s, dirs, dirs)          # 2nd order: the expensive expansion
    .pred <- get("rx_pred_", .s)
    .Dn <- function(.e, .v) symengine::D(.e, symengine::S(.v))
    .sn1 <- function(.j, ...) symengine::S(paste0("rx__sens_", .j, "_BY_", paste(c(...), collapse = "_BY_"), "__"))
    .toRx <- function(.l) rxode2::rxFromSE(.l)
    .g1 <- function(.ex, .p) { .e <- .Dn(.ex, .p); for (.j in .st) .e <- .e + .Dn(.ex, .j) * .sn1(.j, .p); .e }
    .g2 <- function(.ex, .p, .q) { .gq <- .g1(.ex, .q); .e <- .Dn(.gq, .p)
      for (.k in .st) .e <- .e + .Dn(.gq, .k) * .sn1(.k, .p); for (.j in .st) .e <- .e + .Dn(.ex, .j) * .sn1(.j, .p, .q); .e }
    # f2_i_j == f2_j_i, so only the i<=j triangle is differentiated/compiled (~2x off
    # the dominant model-build cost); the reader mirrors A[,i,j]=A[,j,i].
    .P2 <- expand.grid(i = dirs, j = dirs, stringsAsFactors = FALSE)
    .P2 <- .P2[match(.P2$i, dirs) <= match(.P2$j, dirs), , drop = FALSE]
    .fL1 <- vapply(dirs, function(.p) paste0("f1_", .p, "=", .toRx(.g1(.pred, .p))), character(1))
    .fL2 <- vapply(seq_len(nrow(.P2)), function(.r)
      paste0("f2_", .P2$i[.r], "_", .P2$j[.r], "=", .toRx(.g2(.pred, .P2$i[.r], .P2$j[.r]))), character(1))
    .baseOde <- vapply(.st, function(.x) paste0("d/dt(", .x, ")=", .toRx(get(paste0("rx__d_dt_", .x, "__"), .s))), character(1))
    .modTxt <- paste(c(.baseOde, .s1, .s2, paste0("predf=", .toRx(.pred)), .fL1, .fL2), collapse = "\n")
    .modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", .modTxt); .modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", .modTxt)
    list(augMod = rxode2::rxode2(.modTxt), dirs = dirs, ndir = length(dirs), st = .st, P2 = .P2)
  }, error = function(e) NULL)
}

#' Solve the direction-set 2nd-order model for one subject and recover the
#' 3rd-order tensor `Ath` by Shi (2021) central differences of the analytic
#' 2nd-order sensitivities `A` (C++ `shi21CentralWrap`, perturbing one full-param
#' coordinate per direction), then symmetrize.  Returns `list(f, a, A, Ath)` or
#' `NULL` on failure.  (`rxode2::rxExpandSens3_` would give `Ath` analytically but
#' at O(ndir^3) augmented-model compile cost; the Shi FD keeps it at O(ndir^2).)
#' @noRd
.foceiAnalyticSolveSubjectFD3 <- function(aug, params, ev, times, .fdEps = 7e-7, tol = 1e-10) {
  dirs <- aug$dirs; nd <- length(dirs)
  # base solve (f, a, A) via the shared helper; this tier adds the 3rd-order Ath by
  # Shi-differencing A.
  E0 <- .foceiAnalyticSolveFA(aug, params, ev, times, tol = tol); if (is.null(E0)) return(NULL)
  nobs <- length(E0$f); f0 <- as.vector(E0$A)
  # shi21CentralWrap differences the closure at a perturbed full-param vector; it
  # strips names, so re-attach names(params) before mapping back through the solve.
  Aflat <- function(.tt) { E <- .foceiAnalyticSolveFA(aug, setNames(.tt, names(params)), ev, times, tol = tol)
    if (is.null(E) || !all(is.finite(E$A))) return(NULL); as.vector(E$A) }
  Ath <- array(0, c(nobs, nd, nd, nd))
  for (d in seq_len(nd)) {
    idx <- match(dirs[d], names(params))                # this direction's coordinate in params
    if (is.na(idx)) return(NULL)
    sc <- shi21CentralWrap(Aflat, params, f0, idx, .fdEps)  # C++ shi21Central
    if (is.null(sc$gr) || !all(is.finite(sc$gr))) return(NULL)
    Ath[, , , d] <- array(sc$gr, c(nobs, nd, nd))
  }
  # symmetrize over the 3 tensor axes (obs axis 1 held): average the 6 aperm
  # permutations of the tensor axes.  aperm uses the INVERSE-permutation convention,
  # so these are order() of the six S3 permutations -- summed in the same order as the
  # old triple loop (bit-identical).
  pr <- list(c(1,2,3,4), c(1,2,4,3), c(1,3,2,4), c(1,4,2,3), c(1,3,4,2), c(1,4,3,2))
  S <- Reduce(`+`, lapply(pr, function(p) aperm(Ath, p))) / 6
  .out <- list(f = E0$f, a = E0$a, A = E0$A, Ath = S)
  if (!all(is.finite(.out$f)) || !all(is.finite(.out$a)) || !all(is.finite(.out$A)) || !all(is.finite(.out$Ath))) return(NULL)
  .out
}

#' Per-subject observed-information R over structural theta, sigma and Omega, from
#' already-evaluated sensitivities `E$a/A/Ath` + the error model `ef` + the Omega
#' derivatives `omd`.  The inner Hessian / sigma / Omega machinery is eta-indexed
#' (`ei`); the sensitivity slots (N, Tn, dHtD, d2HtDD, and the theta accessors) use
#' the direction index (`di`), and each structural theta differentiates in its own
#' direction-slot via `dirTh` (a mu-ref theta reuses its eta's direction, so a
#' fully mu-referenced model has ndir == neta).  Param order: `nth` theta, `nsg`
#' sigma, then Omega.
#' @noRd
.foceiAnalyticSubjectR <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, omd,
                                   ndir = neta, dirTh = seq_len(nth), Oi = solve(Om),
                                   interaction = 1L) {
  if (identical(as.integer(interaction), 0L))
    return(.foceiAnalyticSubjectRfoce(E, ehat, Om, ef, neta, nth, nsg, sgVar, omd,
                                      ndir = ndir, dirTh = dirTh, Oi = Oi))
  tr <- function(M) sum(diag(M))
  a <- E$a; A <- E$A; Ath <- E$Ath; f <- E$f; y <- E$y
  evf <- function(e) ef$ev(e, f, y)
  rd <- list(r1 = evf(ef$sc$r1), r2 = evf(ef$sc$r2), r3 = evf(ef$sc$r3))
  pf <- list(p = evf(ef$sc$p), p1 = evf(ef$sc$p1), p2 = evf(ef$sc$p2))
  np <- nth + nsg + omd$nom
  ei <- seq_len(neta); di <- seq_len(ndir)
  ae <- a[, ei, drop = FALSE]                          # eta-cols for sigma/Omega
  .dirOf <- function(p) dirTh[p]
  H <- Oi; for (l in ei) for (m in ei) H[l, m] <- H[l, m] + sum(rd$r2 * a[, l] * a[, m] + rd$r1 * A[, l, m])
  N <- matrix(0, neta, ndir); for (l in ei) for (d in di) N[l, d] <- sum(rd$r2 * a[, l] * a[, d] + rd$r1 * A[, l, d])
  HiM <- solve(H)
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(pf$p * a[, l] * a[, m]); Hti <- solve(Ht)
  ouAA <- function(v) { M <- matrix(0, neta, neta); for (l in ei) for (m in ei) M[l, m] <- sum(v * a[, l] * a[, m]); M }
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(pf$p1 * a[, s] * a[, l] * a[, m] + pf$p * A[, l, s] * a[, m] + pf$p * a[, l] * A[, m, s]); D })
  d2HtDD <- lapply(di, function(s) lapply(di, function(t) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(pf$p2 * a[, s] * a[, t] * a[, l] * a[, m] +
      pf$p1 * (A[, s, t] * a[, l] * a[, m] + a[, s] * A[, l, t] * a[, m] + a[, s] * a[, l] * A[, m, t] + a[, t] * A[, l, s] * a[, m] + a[, t] * a[, l] * A[, m, s]) +
      pf$p * (Ath[, l, s, t] * a[, m] + A[, l, s] * A[, m, t] + A[, l, t] * A[, m, s] + a[, l] * Ath[, m, s, t])); D }))
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtD[[l]]), numeric(1))
  Cee <- matrix(0, neta, neta); for (s in ei) for (t in ei)
    Cee[s, t] <- 0.5 * (tr(Hti %*% d2HtDD[[s]][[t]]) - tr(Hti %*% dHtD[[s]] %*% Hti %*% dHtD[[t]]))
  Tn <- array(0, c(neta, ndir, ndir)); for (l in ei) for (s in di) for (t in di)
    Tn[l, s, t] <- sum(rd$r3 * a[, l] * a[, s] * a[, t] + rd$r2 * (A[, l, s] * a[, t] + A[, l, t] * a[, s] + A[, s, t] * a[, l]) + rd$r1 * Ath[, l, s, t])
  Ndd <- function(da, db) sum(rd$r2 * a[, da] * a[, db] + rd$r1 * A[, da, db])
  typ <- function(p) if (p <= nth) "th" else if (p <= nth + nsg) "sg" else "om"
  sgi <- function(p) sgVar[p - nth]
  omc <- function(p) p - nth - nsg
  PVper <- function(p) lapply(ef$per[[sgi(p)]], evf)
  PVpair <- function(aa, bb) { s1 <- sgi(aa); s2 <- sgi(bb); key <- if (paste0(s1, s2) %in% names(ef$pair)) paste0(s1, s2) else paste0(s2, s1)
    lapply(ef$pair[[key]], evf) }
  # Omega enters only the prior (1/2 eta' Omega^-1 eta + 1/2 ln|Omega|) and H~'s
  # +Omega^-1 term, so every Omega-block quantity is an E-basis contraction from
  # `omd` (non-Cholesky variance-covariance derivatives) -- diagonal or block.
  Mcol <- function(p) { t <- typ(p)
    if (t == "th") return(N[, .dirOf(p)]); if (t == "sg") return(as.numeric(crossprod(ae, PVper(p)$rf)))
    as.numeric(omd$dOi[[omc(p)]] %*% ehat) }
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") return(dHtD[[.dirOf(p)]]); if (t == "sg") return(ouAA(PVper(p)$ps))
    omd$dOi[[omc(p)]] }
  d2HtEtaP <- function(p, l) { t <- typ(p)
    if (t == "th") return(d2HtDD[[.dirOf(p)]][[l]]); if (t == "om") return(matrix(0, neta, neta))
    P <- PVper(p); D <- matrix(0, neta, neta); for (s in ei) for (m in ei)
      D[s, m] <- sum(P$pf * a[, l] * a[, s] * a[, m] + P$ps * (A[, s, l] * a[, m] + a[, s] * A[, m, l])); D }
  d2Ht_pp <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "th" && tb == "th") return(d2HtDD[[.dirOf(aa)]][[.dirOf(bb)]])
    if (ta == "om" && tb == "om") return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    if (ta == "om" || tb == "om") return(matrix(0, neta, neta))
    if (ta == "sg" && tb == "sg") return(ouAA(PVpair(aa, bb)$pss))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; d2HtEtaP(sg, .dirOf(thp)) }
  Smat <- function(p) { t <- typ(p)
    if (t == "th") { M <- matrix(0, neta, ndir); for (l in ei) for (s in di) M[l, s] <- Tn[l, .dirOf(p), s]; return(M) }
    if (t == "om") return(omd$dOi[[omc(p)]])
    P <- PVper(p); M <- matrix(0, neta, ndir); for (l in ei) for (s in di) M[l, s] <- sum(P$rff * a[, s] * a[, l] + P$rf * A[, l, s]); M }
  Svec <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "th" && tb == "th") return(Tn[, .dirOf(aa), .dirOf(bb)])
    if (ta == "om" && tb == "om") return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    if (ta == "om" || tb == "om") return(rep(0, neta))
    if (ta == "sg" && tb == "sg") return(as.numeric(crossprod(ae, PVpair(aa, bb)$rfss)))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; Smat(sg)[, .dirOf(thp)] }
  d2Phi <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "th" && tb == "th") return(Ndd(.dirOf(aa), .dirOf(bb)))
    if (ta == "om" && tb == "om") return(0.5 * as.numeric(t(ehat) %*% omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat) + 0.5 * omd$d2LD[omc(aa), omc(bb)])
    if (ta == "om" || tb == "om") return(0)
    if (ta == "sg" && tb == "sg") return(sum(PVpair(aa, bb)$rss))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; as.numeric(crossprod(a[, .dirOf(thp)], PVper(sg)$rf)) }
  npE <- np                                        # full theta + sigma + Omega block
  .Mcols <- lapply(1:npE, Mcol)                     # Mcol(p) depends on ONE index -- cache
  etaP <- matrix(vapply(1:npE, function(p) as.numeric(-HiM %*% .Mcols[[p]]), numeric(neta)), nrow = neta)  # neta x npE (neta==1 safe)
  eta2 <- function(aa, bb) { b <- Svec(aa, bb) + Smat(aa)[, ei, drop = FALSE] %*% etaP[, bb] + Smat(bb)[, ei, drop = FALSE] %*% etaP[, aa]
    for (l in ei) b[l] <- b[l] + as.numeric(t(etaP[, aa]) %*% Tn[l, ei, ei] %*% etaP[, bb]); as.numeric(-HiM %*% b) }
  Cpe <- function(p, l) 0.5 * (tr(Hti %*% d2HtEtaP(p, l)) - tr(Hti %*% dHt_p(p) %*% Hti %*% dHtD[[l]]))
  Cpp <- function(aa, bb) 0.5 * (tr(Hti %*% d2Ht_pp(aa, bb)) - tr(Hti %*% dHt_p(aa) %*% Hti %*% dHt_p(bb)))
  .CpeRow <- lapply(1:npE, function(p) vapply(ei, function(l) Cpe(p, l), numeric(1)))  # length-neta row per p
  R <- matrix(0, npE, npE)
  for (aa in 1:npE) for (bb in aa:npE) {                # R is symmetric -- fill upper, mirror
    dat <- d2Phi(aa, bb) - as.numeric(t(.Mcols[[aa]]) %*% HiM %*% .Mcols[[bb]])
    ld <- Cpp(aa, bb) + sum(.CpeRow[[aa]] * etaP[, bb]) + sum(.CpeRow[[bb]] * etaP[, aa]) +
          as.numeric(t(etaP[, aa]) %*% Cee %*% etaP[, bb]) + sum(Cen * eta2(aa, bb))
    R[aa, bb] <- R[bb, aa] <- dat + ld
  }
  R
}

#' Per-subject observed-information R for FOCE (interaction=0).
#'
#' Unlike FOCEI, the FOCE EBE eta-hat does NOT stationarize the full Laplace
#' objective Phi (the inner problem drops the dR/deta interaction term), so the
#' envelope/Schur shortcut used in [.foceiAnalyticSubjectR] fails.  This uses the
#' GENERAL total-derivative Hessian of F = Phi_full + 0.5 log|H~_FOCE|:
#'
#'   R_ab = F_ab + F_aeta eta_b + F_beta eta_a + eta_a' F_etaeta eta_b + F_eta eta_ab
#'
#' with eta_b, eta_ab the 1st/2nd parameter-derivatives of the EBE obtained by
#' implicit differentiation of the FOCE inner stationarity S_FOCE = sum(q a) +
#' Omega^-1 eta = 0 (q = -eps/R, the interaction-free gradient coefficient; its
#' Jacobian is the FOCE inner Hessian Hf = sum(q' a a' + q A) + Omega^-1 and its
#' 3-tensor Tnf).  Phi keeps the FULL rho derivatives (r1/r2 = Phi_eta / Phi_etaeta
#' and Phi_eta is nonzero here); the log-determinant uses p = 1/R (no interaction).
#' Splitting F into the Phi (data) and 0.5 log|H~| (log-det) parts, the data term is
#' the general (non-envelope) form and the log-det term is the same moving-mode
#' assembly as FOCEI but built from p=1/R.  Same parameter order as
#' [.foceiAnalyticSubjectR].
#' @noRd
.foceiAnalyticSubjectRfoce <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, omd,
                                       ndir = neta, dirTh = seq_len(nth), Oi = solve(Om)) {
  tr <- function(M) sum(diag(M))
  a <- E$a; A <- E$A; Ath <- E$Ath; f <- E$f; y <- E$y
  evf <- function(e) ef$ev(e, f, y)
  rd <- list(r1 = evf(ef$sc$r1), r2 = evf(ef$sc$r2), r3 = evf(ef$sc$r3))   # full rho (Phi)
  qd <- list(q0 = evf(ef$sc$q0), q1 = evf(ef$sc$q1), q2 = evf(ef$sc$q2))   # FOCE inner gradient coef
  pf <- list(p = evf(ef$sc$pF), p1 = evf(ef$sc$pF1), p2 = evf(ef$sc$pF2))  # FOCE determinant p = 1/R
  np <- nth + nsg + omd$nom
  ei <- seq_len(neta); di <- seq_len(ndir); ae <- a[, ei, drop = FALSE]
  .dirOf <- function(p) dirTh[p]

  # ---- Phi (data) pieces: FULL rho derivatives ----
  H <- Oi; for (l in ei) for (m in ei) H[l, m] <- H[l, m] + sum(rd$r2 * a[, l] * a[, m] + rd$r1 * A[, l, m])  # Phi_etaeta
  Ndat <- matrix(0, neta, ndir); for (l in ei) for (d in di) Ndat[l, d] <- sum(rd$r2 * a[, l] * a[, d] + rd$r1 * A[, l, d])  # Phi_(eta,theta)
  gPhi <- as.numeric(Oi %*% ehat); for (l in ei) gPhi[l] <- gPhi[l] + sum(rd$r1 * a[, l])  # Phi_eta (nonzero at eta-hat_FOCE)

  # ---- FOCE inner (EBE) pieces: q-based Jacobian S_eta = Hf and its 3-tensor ----
  Hf <- Oi; for (l in ei) for (m in ei) Hf[l, m] <- Hf[l, m] + sum(qd$q1 * a[, l] * a[, m] + qd$q0 * A[, l, m])  # S_eta = Hf
  HfInv <- solve(Hf)
  Nf <- matrix(0, neta, ndir); for (l in ei) for (d in di) Nf[l, d] <- sum(qd$q1 * a[, l] * a[, d] + qd$q0 * A[, l, d])  # S_(eta,theta)
  Tnf <- array(0, c(neta, ndir, ndir)); for (l in ei) for (s in di) for (t in di)
    Tnf[l, s, t] <- sum(qd$q2 * a[, l] * a[, s] * a[, t] + qd$q1 * (A[, l, s] * a[, t] + A[, l, t] * a[, s] + A[, s, t] * a[, l]) + qd$q0 * Ath[, l, s, t])  # S_etaeta

  # ---- determinant H~ (p = 1/R) pieces (same moving-mode assembly as FOCEI) ----
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(pf$p * a[, l] * a[, m]); Hti <- solve(Ht)
  ouAA <- function(v) { M <- matrix(0, neta, neta); for (l in ei) for (m in ei) M[l, m] <- sum(v * a[, l] * a[, m]); M }
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(pf$p1 * a[, s] * a[, l] * a[, m] + pf$p * A[, l, s] * a[, m] + pf$p * a[, l] * A[, m, s]); D })
  d2HtDD <- lapply(di, function(s) lapply(di, function(t) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(pf$p2 * a[, s] * a[, t] * a[, l] * a[, m] +
      pf$p1 * (A[, s, t] * a[, l] * a[, m] + a[, s] * A[, l, t] * a[, m] + a[, s] * a[, l] * A[, m, t] + a[, t] * A[, l, s] * a[, m] + a[, t] * a[, l] * A[, m, s]) +
      pf$p * (Ath[, l, s, t] * a[, m] + A[, l, s] * A[, m, t] + A[, l, t] * A[, m, s] + a[, l] * Ath[, m, s, t])); D }))
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtD[[l]]), numeric(1))
  Cee <- matrix(0, neta, neta); for (s in ei) for (t in ei)
    Cee[s, t] <- 0.5 * (tr(Hti %*% d2HtDD[[s]][[t]]) - tr(Hti %*% dHtD[[s]] %*% Hti %*% dHtD[[t]]))

  typ <- function(p) if (p <= nth) "th" else if (p <= nth + nsg) "sg" else "om"
  sgi <- function(p) sgVar[p - nth]; omc <- function(p) p - nth - nsg
  perR <- function(p) lapply(ef$per[[sgi(p)]][c("rf", "rff")], evf)                 # rho sigma partials (Phi)
  perQ <- function(p) list(qs = evf(ef$per[[sgi(p)]]$qs), qsf = evf(ef$per[[sgi(p)]]$qsf))  # FOCE inner sigma
  perP <- function(p) list(ps = evf(ef$per[[sgi(p)]]$psF), pf = evf(ef$per[[sgi(p)]]$pfF))   # FOCE det sigma
  pairKey <- function(aa, bb) { s1 <- sgi(aa); s2 <- sgi(bb); if (paste0(s1, s2) %in% names(ef$pair)) paste0(s1, s2) else paste0(s2, s1) }
  pairR <- function(aa, bb) lapply(ef$pair[[pairKey(aa, bb)]][c("rss", "rfss")], evf)
  pairQ <- function(aa, bb) list(qss = evf(ef$pair[[pairKey(aa, bb)]]$qss))
  pairP <- function(aa, bb) list(pss = evf(ef$pair[[pairKey(aa, bb)]]$pssF))

  # ---- data-side accessors (Phi partials) ----
  McolData <- function(p) { t <- typ(p)                                            # Phi_(eta,p)
    if (t == "th") return(Ndat[, .dirOf(p)]); if (t == "sg") return(as.numeric(crossprod(ae, perR(p)$rf)))
    as.numeric(omd$dOi[[omc(p)]] %*% ehat) }
  d2Phi <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)                          # Phi_(p,p) at fixed eta
    if (ta == "th" && tb == "th") return(sum(rd$r2 * a[, .dirOf(aa)] * a[, .dirOf(bb)] + rd$r1 * A[, .dirOf(aa), .dirOf(bb)]))
    if (ta == "om" && tb == "om") return(0.5 * as.numeric(t(ehat) %*% omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat) + 0.5 * omd$d2LD[omc(aa), omc(bb)])
    if (ta == "om" || tb == "om") return(0)
    if (ta == "sg" && tb == "sg") return(sum(pairR(aa, bb)$rss))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; as.numeric(crossprod(a[, .dirOf(thp)], perR(sg)$rf)) }

  # ---- EBE-side accessors (S_FOCE partials) ----
  McolEBE <- function(p) { t <- typ(p)                                             # S_p
    if (t == "th") return(Nf[, .dirOf(p)]); if (t == "sg") return(as.numeric(crossprod(ae, perQ(p)$qs)))
    as.numeric(omd$dOi[[omc(p)]] %*% ehat) }
  SmatEBE <- function(p) { t <- typ(p)                                             # S_(p,eta)
    if (t == "th") { M <- matrix(0, neta, ndir); for (l in ei) for (s in di) M[l, s] <- Tnf[l, .dirOf(p), s]; return(M) }
    if (t == "om") return(omd$dOi[[omc(p)]])
    P <- perQ(p); M <- matrix(0, neta, ndir); for (l in ei) for (s in di) M[l, s] <- sum(P$qsf * a[, s] * a[, l] + P$qs * A[, l, s]); M }
  SvecEBE <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)                       # S_(p,p)
    if (ta == "th" && tb == "th") return(Tnf[, .dirOf(aa), .dirOf(bb)])
    if (ta == "om" && tb == "om") return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    if (ta == "om" || tb == "om") return(rep(0, neta))
    if (ta == "sg" && tb == "sg") return(as.numeric(crossprod(ae, pairQ(aa, bb)$qss)))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; SmatEBE(sg)[, .dirOf(thp)] }

  # ---- determinant-side accessors (H~ partials, p = 1/R) ----
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") return(dHtD[[.dirOf(p)]]); if (t == "sg") return(ouAA(perP(p)$ps)); omd$dOi[[omc(p)]] }
  d2HtEtaP <- function(p, l) { t <- typ(p)
    if (t == "th") return(d2HtDD[[.dirOf(p)]][[l]]); if (t == "om") return(matrix(0, neta, neta))
    P <- perP(p); D <- matrix(0, neta, neta); for (s in ei) for (m in ei)
      D[s, m] <- sum(P$pf * a[, l] * a[, s] * a[, m] + P$ps * (A[, s, l] * a[, m] + a[, s] * A[, m, l])); D }
  d2Ht_pp <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "th" && tb == "th") return(d2HtDD[[.dirOf(aa)]][[.dirOf(bb)]])
    if (ta == "om" && tb == "om") return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    if (ta == "om" || tb == "om") return(matrix(0, neta, neta))
    if (ta == "sg" && tb == "sg") return(ouAA(pairP(aa, bb)$pss))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; d2HtEtaP(sg, .dirOf(thp)) }
  Cpe <- function(p, l) 0.5 * (tr(Hti %*% d2HtEtaP(p, l)) - tr(Hti %*% dHt_p(p) %*% Hti %*% dHtD[[l]]))
  Cpp <- function(aa, bb) 0.5 * (tr(Hti %*% d2Ht_pp(aa, bb)) - tr(Hti %*% dHt_p(aa) %*% Hti %*% dHt_p(bb)))

  npE <- np
  .McD <- lapply(1:npE, McolData)                                                  # Phi_(eta,p), cached
  etaP <- matrix(vapply(1:npE, function(p) as.numeric(-HfInv %*% McolEBE(p)), numeric(neta)), nrow = neta)  # eta_p = -Hf^-1 S_p
  eta2 <- function(aa, bb) {                                                        # eta_ab (2nd EBE deriv)
    b <- SvecEBE(aa, bb) + SmatEBE(aa)[, ei, drop = FALSE] %*% etaP[, bb] + SmatEBE(bb)[, ei, drop = FALSE] %*% etaP[, aa]
    for (l in ei) b[l] <- b[l] + as.numeric(t(etaP[, aa]) %*% Tnf[l, ei, ei] %*% etaP[, bb]); as.numeric(-HfInv %*% b) }
  .CpeRow <- lapply(1:npE, function(p) vapply(ei, function(l) Cpe(p, l), numeric(1)))
  R <- matrix(0, npE, npE)
  for (aa in 1:npE) for (bb in aa:npE) {                # R is symmetric -- fill upper, mirror
    e_ab <- eta2(aa, bb)
    # data term (general, non-envelope): F_ab + F_aeta eta_b + F_beta eta_a + eta_a' F_etaeta eta_b + F_eta eta_ab, Phi part
    dat <- d2Phi(aa, bb) + sum(.McD[[aa]] * etaP[, bb]) + sum(.McD[[bb]] * etaP[, aa]) +
           as.numeric(t(etaP[, aa]) %*% H %*% etaP[, bb]) + sum(gPhi * e_ab)
    # log-determinant term (moving mode), 0.5 log|H~| part
    ld <- Cpp(aa, bb) + sum(.CpeRow[[aa]] * etaP[, bb]) + sum(.CpeRow[[bb]] * etaP[, aa]) +
          as.numeric(t(etaP[, aa]) %*% Cee %*% etaP[, bb]) + sum(Cen * e_ab)
    R[aa, bb] <- R[bb, aa] <- dat + ld
  }
  R
}

#' Base subject solve: `f` plus the 1st/2nd analytic sensitivities (`a`, `A`), no
#' 3rd-order Shi tensor.  Shared by [.foceiAnalyticSolveSubjectFD3] (which adds `Ath`)
#' and the FOCE EBE re-solve (which needs only `a`/`A`).  Muffles benign solver
#' warnings (a real error returns `NULL` -> FD fallback); the nrow guard bails when an
#' EVID==2/covariate-update row shares an obs timestamp (would misalign f against y).
#' @noRd
.foceiAnalyticSolveFA <- function(aug, params, ev, times, tol = 1e-10) {
  dirs <- aug$dirs; nd <- length(dirs)
  .d <- tryCatch(withCallingHandlers(
      as.data.frame(rxode2::rxSolve(aug$augMod, params = params, ev,
          returnType = "data.frame", atol = tol, rtol = tol)),
      warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL)
  if (is.null(.d)) return(NULL); .d <- .d[.d$time %in% times, , drop = FALSE]
  if (nrow(.d) == 0L || nrow(.d) != length(times) ||
        !all(c("predf", paste0("f1_", dirs)) %in% names(.d))) return(NULL)
  no <- nrow(.d)
  a <- matrix(vapply(dirs, function(q) .d[[paste0("f1_", q)]], numeric(no)), no, nd)
  A <- array(0, c(no, nd, nd))
  for (r in seq_len(nrow(aug$P2))) {                   # only i<=j emitted -> mirror to i>j
    .ii <- match(aug$P2$i[r], dirs); .jj <- match(aug$P2$j[r], dirs)
    .v2 <- .d[[paste0("f2_", aug$P2$i[r], "_", aug$P2$j[r])]]
    A[, .ii, .jj] <- .v2; A[, .jj, .ii] <- .v2
  }
  list(f = .d$predf, a = a, A = A)
}

#' Re-solve one subject's EBE to the FOCE inner stationary point S_FOCE = sum(q a)
#' + Omega^-1 eta = 0 (q = -eps/R) via Newton on the FOCE inner Hessian
#' Hf = sum(q' a a' + q A) + Omega^-1, starting from the stored eta `eta0`.
#' nlmixr's stored FOCE-combined EBEs do NOT satisfy S_FOCE=0 (an estimation-side
#' inconsistency), so R must be formed at the re-solved eta.  For additive/FOCEI the
#' stored eta is already stationary (|S_FOCE| < `skip`) -> returns `eta0` unchanged
#' (byte no-op).  `NULL` on a solve/Newton failure -> caller falls back to FD.
#' @noRd
.foceiAnalyticFoceEbe <- function(aug, th, eta0, s, times, y, etav, ef, Oi, neta, tol,
                                  maxit = 30L, skip = 1e-3, conv = 1e-9) {
  ei <- seq_len(neta)
  .SH <- function(eta) {                               # FOCE S_FOCE and its Jacobian Hf at eta
    E <- .foceiAnalyticSolveFA(aug, c(th, setNames(eta, etav)), s, times, tol = tol)
    if (is.null(E)) return(NULL)
    q0 <- ef$ev(ef$sc$q0, E$f, y); q1 <- ef$ev(ef$sc$q1, E$f, y)
    S <- as.numeric(Oi %*% eta); for (l in ei) S[l] <- S[l] + sum(q0 * E$a[, l])
    Hf <- Oi; for (l in ei) for (m in ei) Hf[l, m] <- Hf[l, m] + sum(q1 * E$a[, l] * E$a[, m] + q0 * E$A[, l, m])
    list(S = S, Hf = Hf)
  }
  eta <- eta0
  sh <- .SH(eta); if (is.null(sh)) return(NULL)
  if (max(abs(sh$S)) < skip) return(eta0)              # already FOCE-stationary (additive/FOCEI) -> no-op
  for (it in seq_len(maxit)) {
    step <- tryCatch(solve(sh$Hf, sh$S), error = function(e) NULL)
    if (is.null(step)) return(NULL)
    eta <- eta - step
    sh <- .SH(eta); if (is.null(sh)) return(NULL)
    if (max(abs(sh$S)) < conv) break
  }
  if (max(abs(sh$S)) >= conv) return(NULL)               # Newton did not converge -> FD fallback
  eta
}

#' Full analytic FOCEI covariance (theta + sigma + Omega) for a fitted object, or
#' `NULL` when out of scope / the augmented solve fails.  Standalone test oracle.
#' @param fit a fitted nlmixr2 focei object
#' @return list(cov, se, R, params, method) or `NULL`
#' @noRd
foceiCovAnalytic <- function(fit) {
  ui <- fit$finalUi
  if (!.hasRxSens()) return(NULL)            # need rxExpandSens2_ + symengine
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(NULL)             # FO/FOI
  interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))                   # 1 FOCEI / 0 FOCE
  if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L) return(NULL)                  # AGQ
  # IOV is out of scope; derive the flag from THIS fit, not the process-global
  # .uiIovEnv (which reflects the LAST-preprocessed model).  An occasion (IOV) eta
  # carries a non-"id" `condition` -- the same signal .uiApplyIov keys on.
  .idf0 <- ui$iniDf
  if (any(!is.na(.idf0$condition) & .idf0$condition != "id" & is.na(.idf0$err))) return(NULL)  # IOV
  if (!is.null(fit$dataSav$CENS) && any(fit$dataSav$CENS != 0, na.rm = TRUE)) return(NULL) # censored
  ef <- .foceiAnalyticErrFull(ui)
  if (is.null(ef)) return(NULL)                     # unsupported error model -> fallback

  ini <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)                    # theta <-> eta pairing
  etaNames <- .map$etaNames
  neta <- length(etaNames)
  if (neta == 0L) return(NULL)
  thetaForEta <- .map$thetaForEta
  if (anyNA(thetaForEta)) return(NULL)             # non-mu-ref eta -> out of scope
  if (any(.iniIsFixed(ini, thetaForEta))) return(NULL)  # fixed structural theta breaks eta indexing
  keep <- !.iniIsFixed(ini, ef$sgName); ef$sgVar <- ef$sgVar[keep]; ef$sgName <- ef$sgName[keep]  # drop fixed sigma
  Om <- fit$omega
  pairs <- .foceiOmegaPairs(Om, ini)               # free Omega lower-triangle (declared blocks)
  omd <- .omegaVarCovDeriv(Om, pairs)

  # Uniform direction assembly (shared with the production covType="analytic" hook).
  .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta)
  if (is.null(.dir)) return(NULL)
  thStruct <- .dir$thStruct; dirs <- .dir$dirs; dirTh <- .dir$dirTh
  ndir <- .dir$ndir; nth <- .dir$nth

  thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]
  thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
  # converged estimates (fit$theta), not ui$iniDf$est which can hold initials
  .thv <- fit$theta[thRows$name]
  if (anyNA(.thv)) .thv <- thRows$est
  th <- setNames(as.numeric(.thv), paste0("THETA_", seq_len(nrow(thRows)), "_"))
  ebes <- as.matrix(fit$eta[, etaNames, drop = FALSE])
  nsg <- length(ef$sgVar)

  R <- .foceiAnalyticAssembleR(ui, th, ebes, fit$eta$ID, fit$dataSav, Om, ef, neta, nth, nsg, omd,
                               dirs = dirs, dirTh = dirTh, ndir = ndir,
                               solveTol = .foceiAnalyticSolveTol(ui), interaction = interaction)
  if (is.null(R)) return(NULL)
  cov <- tryCatch(solve(R), error = function(e) NULL)
  if (is.null(cov)) return(NULL)
  nm <- c(thStruct, ef$sgName, .foceiOmegaCovNames(pairs, thetaForEta))
  dimnames(R) <- dimnames(cov) <- list(nm, nm)
  list(cov = cov, se = setNames(suppressWarnings(sqrt(diag(cov))), nm),  # NaN flags non-PD
       R = R, params = nm, method = "analytic")
}
