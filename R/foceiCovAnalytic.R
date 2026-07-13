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
  .ret$covMethod <- "analytic"           # report the analytic observed information (not "r")
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

#' Emit the reason the analytic (covType="analytic") observed-information R-matrix
#' is unavailable for this model and return `NULL` so the caller drops to the
#' finite-difference Hessian.  Only reached in opted-in paths -- the native
#' `.foceiCalcRanalytic` hook checks `covType=="analytic"` first and
#' `foceiCovAnalytic` is an explicit call -- so nothing prints for the default
#' `covType="fd"`.  Kept as a plain `message()` (informational, not a warning):
#' falling back to FD is a valid result, not an error.
#' @param reason short human-readable phrase naming the out-of-scope feature
#' @return `NULL`, invisibly usable as `return(.foceiAnalyticFallback(...))`
#' @noRd
.foceiAnalyticFallback <- function(reason) {
  message("covType=\"analytic\": ", reason,
          " is out of analytic-covariance scope; using the finite-difference covariance instead")
  NULL
}

#' Transform the observed DV onto the rx_pred_ (transformed) scale for a both-sides
#' transform, using the solved per-observation transform parameters `trans`
#' (lambda/yj/low/hi, each length nobs).  Runs the SAME C++ transform the inner problem
#' uses (`_nlmixr2est_powerD` == `tbs()`), which requires all arguments the same length --
#' the per-observation vectors already are, so a multi-endpoint mix of transforms is
#' handled directly.  Returns `dv` unchanged when `trans` is NULL.
#' @noRd
.foceiAnalyticTbsY <- function(dv, trans) {
  if (is.null(trans)) return(dv)
  .Call(`_nlmixr2est_powerD`, as.double(dv), as.double(trans$lambda),
        as.integer(trans$yj), as.double(trans$low), as.double(trans$hi))
}

#' Per-observation transform value lambda-derivative dy'/dlambda of the DV
#' (`d tbs(DV,lambda)/dlambda`), used to build the residual DV-transform chain
#' when lambda is estimated.  Same per-obs `trans` contract as `.foceiAnalyticTbsY`.
#' @noRd
.foceiAnalyticDvSensLambda <- function(dv, trans) {
  if (is.null(trans)) return(rep(0, length(dv)))
  .Call(`_nlmixr2est_powerDLambda`, as.double(dv), as.double(trans$lambda),
        as.integer(trans$yj), as.double(trans$low), as.double(trans$hi))
}

#' Per-observation d2 tbs(DV,lambda)/dlambda2 (DV second lambda-derivative).
#' @noRd
.foceiAnalyticDvSensLambda2 <- function(dv, trans) {
  if (is.null(trans)) return(rep(0, length(dv)))
  .Call(`_nlmixr2est_powerDLambda2`, as.double(dv), as.double(trans$lambda),
        as.integer(trans$yj), as.double(trans$low), as.double(trans$hi))
}

#' Per-observation d log|dy'/dDV| / dlambda -- the transform Jacobian's
#' lambda-derivative, the extra term the OFV gradient's lambda column carries
#' (-2 * sum over obs); it cancels in the observed-information covariance.
#' @noRd
.foceiAnalyticJacLambda <- function(dv, trans) {
  if (is.null(trans)) return(rep(0, length(dv)))
  .Call(`_nlmixr2est_powerDL`, as.double(dv), as.double(trans$lambda),
        as.integer(trans$yj), as.double(trans$low), as.double(trans$hi))
}

#' Count literal occurrences of each of `.names` on the model-expression RHS.
#'
#' Shared by the eta-occurrence guard (a random effect reused across parameters breaks
#' the covariate/theta direction reuse) and the structural-theta guard (a structural
#' parameter appearing in more than one mu-expression breaks the eta-less covariate
#' reuse -- `df/dtheta` would then gain a path `df/db` lacks).  Returns a named integer
#' vector; `2L` (i.e. ">1", conservative -> no reuse) when `lstExpr` is unavailable.
#' @noRd
.foceiNameOccurrence <- function(ui, .names) {
  .lst <- tryCatch(ui$lstExpr, error = function(e) NULL)
  .count <- function(.e, .nm) {
    if (is.name(.e)) return(as.integer(identical(as.character(.e), .nm)))
    if (is.call(.e)) return(sum(vapply(as.list(.e)[-1L], .count, integer(1), .nm)))
    0L
  }
  .isAssign <- function(.ex) is.call(.ex) && length(.ex) == 3L &&
    (identical(.ex[[1]], as.name("<-")) || identical(.ex[[1]], as.name("=")))
  stats::setNames(vapply(as.character(.names), function(.n)
    if (is.null(.lst)) 2L else sum(vapply(.lst, function(.ex)
      if (.isAssign(.ex)) .count(.ex[[3]], .n) else 0L, integer(1))), integer(1)),
    as.character(.names))
}
#' Literal occurrence count of each diagonal eta across the model RHS.
#'
#' The sensitivity direction-reuse identities -- a mu-ref theta reusing its eta's
#' direction (`df/dtheta = df/deta`) and a mu-ref covariate coefficient scaling it
#' (`df/db = cov*df/deta`) -- are exact only when the eta enters the model through a
#' single additive position.  An eta shared across parameters (e.g. `eta.cl` in both
#' `cl` and `v`) makes `df/dtheta` (one parameter) differ from `df/deta` (all of
#' them), so the reuse would be wrong.  rxode2's mu-reference frames still map such an
#' eta, so its multiplicity is counted here.  Named by the diagonal-eta names in
#' `neta1` order; `2L` (conservative ">1") when `lstExpr` is unavailable.
#' @noRd
.foceiEtaOccurrence <- function(ui) {
  .idf <- ui$iniDf
  .etaRows <- .idf[!is.na(.idf$neta1) & .idf$neta1 == .idf$neta2, , drop = FALSE]
  .etaRows <- .etaRows[order(.etaRows$neta1), , drop = FALSE]
  .foceiNameOccurrence(ui, as.character(.etaRows$name))
}

#' Scope-relevant structural thetas + sensitivity direction map: one direction per
#' eta (ETA_i_) plus one per non-mu-ref structural theta (THETA_j_); a mu-ref theta
#' reuses its eta's direction (so a fully mu-ref model has ndir == neta).  Sigma and
#' IOV SD thetas (`sgName`/`iovVars`) are not structural.  `sharedEta` is a logical
#' flag per diagonal eta (in `neta1` order) marking an eta reused across parameters,
#' for which the theta reuse is invalid.  `NULL` if none remain.
#' @noRd
.foceiAnalyticDirections <- function(ini, thetaForEta, sgName, neta, iovVars = character(0),
                                     sharedEta = logical(0)) {
  thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]
  thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
  thStructRows <- thRows[!.iniIsFixed(ini, thRows$name) & !(thRows$name %in% sgName) &
                           !(thRows$name %in% iovVars), , drop = FALSE]
  thStruct <- thStructRows$name
  nth <- length(thStruct)
  if (nth == 0L) return(.foceiAnalyticFallback("a model with no estimated structural (fixed-effect) parameter"))
  etaDirs <- paste0("ETA_", seq_len(neta), "_")
  dirTh <- integer(nth); nonMuTheta <- character(0)
  for (p in seq_len(nth)) {
    .mt <- which(thetaForEta == thStruct[p])        # eta(s) this theta mu-references
    if (length(.mt) > 1L)                           # shared by >1 eta -> summing routes NYI, fall back to FD
      return(.foceiAnalyticFallback("a structural parameter mu-referenced by more than one random effect"))
    k <- if (length(.mt) == 1L) .mt else NA_integer_
    # A mu-ref theta reuses its eta's state-sensitivity direction for free (df/dtheta = df/deta,
    # since theta and eta enter the mu identically) -- but ONLY when the eta enters the model in a
    # single additive position.  A shared eta (e.g. eta.cl in both cl and v) makes df/dtheta (this
    # parameter only) differ from df/deta (all parameters the eta appears in), so the reuse is
    # invalid; treat the theta as non-mu and give it its OWN true-sensitivity direction (the eta
    # keeps its own, correct, direction).  Stays analytic -- just one extra direction for that
    # theta, only in the rare shared-eta case.
    .reuse <- !is.na(k) && !(k <= length(sharedEta) && isTRUE(sharedEta[k]))
    if (.reuse) {
      dirTh[p] <- k                                 # reuse eta k's direction (free)
    } else {
      nonMuTheta <- c(nonMuTheta, paste0("THETA_", thStructRows$ntheta[p], "_"))
      dirTh[p] <- neta + length(nonMuTheta)         # its own true-sensitivity direction
    }
  }
  dirs <- c(etaDirs, nonMuTheta)
  # Estimated boxCox/yeoJohnson lambda is an error parameter that (unlike an ordinary
  # sigma with df/dsigma=0) has a NONZERO prediction sensitivity df'/dlambda=rxTBSdL(f,
  # lambda), so it enters as its OWN theta-like direction (appended to the theta block):
  # the kernels then handle its pred + variance sensitivities exactly like a theta, and
  # the DV-transform residual chain (dy'/dlambda) + the -2 log|J| Jacobian are applied on
  # top (in the driver / kernel).  Fixed lambda is filtered out here (static transform).
  .lamRows <- ini[ini$name %in% sgName & !is.na(ini$err) &
                    ini$err %in% c("boxCox", "yeoJohnson") & !.iniIsFixed(ini, ini$name), , drop = FALSE]
  lamNames <- character(0); lamDir <- integer(0)
  if (nrow(.lamRows) > 0L) {
    .lamRows <- .lamRows[order(.lamRows$ntheta), , drop = FALSE]
    for (r in seq_len(nrow(.lamRows))) {
      dirs <- c(dirs, paste0("THETA_", .lamRows$ntheta[r], "_"))
      lamNames <- c(lamNames, .lamRows$name[r]); lamDir <- c(lamDir, length(dirs))
    }
    thStruct <- c(thStruct, lamNames); dirTh <- c(dirTh, lamDir); nth <- nth + length(lamNames)
  }
  lamIdx <- if (length(lamNames)) seq(nth - length(lamNames) + 1L, nth) else integer(0)  # 1-based theta-block positions
  # (f,R) covariance: fold the residual-error parameters (sigmas) into the DIRECTION set.
  # The NONMEM sigma (the unit N(0,1) residual variate) is FIXED to 1, so the estimated
  # error parameters enter ONLY through the weight/variance rx_r_ and NOT through the
  # prediction -- i.e. a sigma direction has df/dsigma=0 (a=A=Ath=0) and only carries the
  # variance sensitivities dR/dsigma (aR/AR/AthR).  Treating each sigma as its own direction
  # therefore lets the SAME 1st/2nd-order solve + Shi-FD 3rd-order machinery supply every
  # sigma derivative (including the 3rd-order mixed sigma-direction terms the determinant
  # block needs) with no special-casing.  NOTE: if the FOCEI methods are ever generalized so
  # the residual sigma itself is estimated (not fixed to 1), that estimated sigma would enter
  # the variance multiplicatively (R -> R*sigma^2) and would need its own chain here.
  .sgRows <- thRows[thRows$name %in% sgName & !(thRows$name %in% lamNames), , drop = FALSE]
  dirSg <- if (nrow(.sgRows) > 0L) seq_len(nrow(.sgRows)) + length(dirs) else integer(0)
  dirsCov <- c(dirs, if (nrow(.sgRows) > 0L) paste0("THETA_", .sgRows$ntheta, "_") else character(0))
  list(thStruct = thStruct, thStructRows = thStructRows, dirs = dirs,
       dirTh = dirTh, ndir = length(dirs), nth = nth,
       dirsCov = dirsCov, ndirCov = length(dirsCov), sgName = .sgRows$name,
       lamNames = lamNames, lamDir = lamDir, lamIdx = lamIdx,
       dirP = c(dirTh, dirSg))                          # every non-Omega param -> a direction
}

#' Sum the per-subject general (f,R) observed-information R over subjects (FOCEI).
#' Builds the augmented model over the covariance direction set (etas + non-mu thetas +
#' SIGMA thetas), Shi-FDs the 3rd-order tensors (Ath + AthR), and assembles via
#' `.foceiAnalyticSubjectRFR`.  Handles any variance structure; `NULL` on failure.
#' @noRd
.foceiAnalyticAssembleRFR <- function(ui, th, ebes, ids, data, Om, ef, neta, ndirP, dirP, omd,
                                      dirsCov, ndirCov, startedEnv = NULL, solveTol = 1e-10,
                                      interaction = 1L, foceType = 0L, lamDir = integer(0)) {
  # Route A (default ON; set FOCEI_NO_RSIG=1 to opt out): build the gradient's `dirs` model
  # (drop the sigma directions) so the (f,R) cov SHARES the gradient's already-compiled model
  # -- one augmented-model compile instead of two, and the sigma directions' wasted (all-zero)
  # state-sensitivity columns vanish.  The sigma tensor slots (aR/AR/AthR) are rebuilt from the
  # model's own rsig outputs in SolveAllFD3 (.foceiAnalyticExpandSigma), which expands the solved
  # E back to ndirCov -- numerically identical to carrying explicit sigma directions.
  .rsigA <- !nzchar(Sys.getenv("FOCEI_NO_RSIG"))
  .erN <- ui$iniDf$ntheta[!is.na(ui$iniDf$err) & !(ui$iniDf$err %in% c("boxCox", "yeoJohnson"))]
  .sigDirs <- if (.rsigA) intersect(dirsCov, paste0("THETA_", .erN, "_")) else character(0)
  .dirsF <- dirsCov[!(dirsCov %in% .sigDirs)]
  # Flag BEFORE building the augmented model: .foceiAnalyticAugModelDirs loads a new
  # compiled model, which calls rxSolveFree() on the fit's global solve.  If the build
  # then declines (is.null(am)/ndir mismatch/no rx_r_), the C++ hook must still restore the
  # freed solve before the FD fallback -- otherwise the FD sandwich solves against freed memory.
  if (!is.null(startedEnv)) assign(".analyticStarted", TRUE, startedEnv)
  am <- .foceiAnalyticAugModelDirs(ui, if (.rsigA) .dirsF else dirsCov)
  if (is.null(am) || am$ndir != (ndirCov - length(.sigDirs)) || !isTRUE(am$hasRvar)) return(NULL)
  np <- ndirP + omd$nom; Oi <- solve(Om)
  etav <- paste0("ETA_", seq_len(neta), "_")
  .foce <- identical(as.integer(interaction), 0L)      # FOCE re-solves EBEs to S_FOCE=0
  .fp <- identical(as.integer(foceType), 1L)           # foce+ keeps the live conditional R
  .byId <- split(data, as.character(data$ID))
  .idCode <- if (is.factor(ids)) as.integer(ids) else match(ids, sort(unique(ids)))
  R <- matrix(0, np, np)
  if (.foce) {
    # FOCE cov: the per-subject eta=0 population solve + EBE re-solve + 3rd-order Shi FD3 stay
    # in R (inherently per-subject), then ONE OpenMP C++ call (foceiRAllFoceFR_) sums the
    # observed information over subjects.  The frozen-R0 sensitivities are resolved per subject
    # (nonmem: aRe/ARe=0, aRc/ARc/R0 from E0; foce+: all from the eta-hat solve E).
    nsub <- length(ids); Elist <- vector("list", nsub); E0list <- vector("list", nsub)
    eta0list <- vector("list", nsub); nobsAll <- integer(nsub)
    # BATCHED (f,R) FOCE/foce+ solves: the eta=0 population solve (nonmem frozen R0) is batched,
    # the EBE re-solve stays per-subject (Newton), then ONE batched SolveAllFD3 delivers f/a/A/Ath
    # AND R/aR/AR for all subjects (withR=FALSE: FOCE reads R/aR/AR from the base solve + Ath, but
    # never AthR).  foce+ (foceType=1) keeps the live R (no eta=0 solve).  Per-subject Shi fallback.
    .obsAll <- lapply(seq_len(nsub), function(i) { .s <- .byId[[as.character(.idCode[i])]]
      if (is.null(.s) || nrow(.s) == 0L) NULL else .s[.s$EVID == 0, , drop = FALSE] })
    if (any(vapply(.obsAll, is.null, logical(1L)))) return(NULL)
    .obsT <- lapply(.obsAll, function(.o) .o$TIME)
    E0all <- if (!.fp) .foceiAnalyticSolveAll(am, th, matrix(0, nsub, neta), .idCode, data, .obsT, solveTol) else NULL
    if (!.fp && is.null(E0all)) return(NULL)
    if (.rsigA && !.fp) {                                # Route A: rebuild frozen-R0 sigma slots (aR/AR) from rsig
      .ns0 <- ncol(E0all[[1L]]$Rsig)
      if (!is.null(.ns0) && .ns0 > 0L) E0all <- lapply(E0all, function(.E0) .foceiAnalyticExpandSigma(.E0, .ns0, neta, NULL, NULL))
    }
    eta0Mat <- .foceiAnalyticFoceEbeBatch(am, th, ebes, .idCode, data, .obsAll, .obsT, etav, Oi, neta,
                                          solveTol, foceType = foceType, E0all = E0all)   # batched EBE re-solve
    if (is.null(eta0Mat)) return(NULL)
    .batch <- !nzchar(Sys.getenv("FOCEI_NO_FD3_BATCH"))
    .EsAll <- if (.batch) .foceiAnalyticSolveAllFD3(am, th, eta0Mat, .idCode, data, .obsT, tol = solveTol, withR = FALSE) else NULL
    if (.batch && is.null(.EsAll)) .batch <- FALSE
    for (i in seq_len(nsub)) {
      s <- .byId[[as.character(.idCode[i])]]; obs <- .obsAll[[i]]; eta0 <- eta0Mat[i, ]
      E0 <- if (.fp) NULL else E0all[[i]]
      E <- if (.batch) .EsAll[[i]] else .foceiAnalyticSolveSubjectFD3(am, c(th, setNames(eta0, etav)), s, obs$TIME, tol = solveTol, withR = FALSE)
      if (is.null(E)) return(NULL)
      if (isTRUE(ef$canVanish)) { .fa <- abs(E$f)
        if (any(!is.finite(.fa)) || min(.fa) < 1e-6 * max(.fa))
          return(.foceiAnalyticFallback("pure proportional error with a near-zero model prediction")) }
      E$y <- .foceiAnalyticTbsY(obs$DV, E$trans)
      # E0list[i] <- list(E0) keeps a NULL slot (foce+) without shrinking the list.
      Elist[[i]] <- E; E0list[i] <- list(E0); eta0list[[i]] <- eta0; nobsAll[i] <- length(E$f)
    }
    .fpG <- identical(as.integer(foceType), 1L) || is.null(E0list[[1L]])
    totObs <- sum(nobsAll); off <- c(0L, cumsum(nobsAll)); nd2 <- ndirCov * ndirCov
    aB <- matrix(0, totObs, ndirCov); aReB <- matrix(0, totObs, ndirCov); aRcB <- matrix(0, totObs, ndirCov)
    AB <- array(0, c(totObs, ndirCov, ndirCov)); AReB <- array(0, c(totObs, ndirCov, ndirCov)); ARcB <- array(0, c(totObs, ndirCov, ndirCov))
    AthB <- array(0, c(totObs, neta, nd2))
    fB <- numeric(totObs); yB <- numeric(totObs); R0B <- numeric(totObs); ehatB <- matrix(0, nsub, neta)
    dvSensB <- if (length(lamDir)) matrix(0, totObs, ndirCov) else matrix(0, totObs, 0L)
    dvSens2B <- dvSensB
    .hasCens <- (!is.null(data$CENS) && any(data$CENS != 0, na.rm = TRUE)) ||
      (!is.null(data$LIMIT) && any(is.finite(data$LIMIT)))
    censB <- if (.hasCens) integer(totObs) else integer(0)
    limB <- if (.hasCens) rep(NA_real_, totObs) else numeric(0)
    for (i in seq_len(nsub)) {
      E <- Elist[[i]]; E0 <- E0list[[i]]; no <- nobsAll[i]; rows <- (off[i] + 1L):off[i + 1L]
      aB[rows, ] <- E$a; AB[rows, , ] <- E$A; AthB[rows, , ] <- array(E$Ath, c(no, neta, nd2))
      fB[rows] <- E$f; yB[rows] <- E$y
      if (.fpG) { R0B[rows] <- E$R; aReB[rows, ] <- E$aR; aRcB[rows, ] <- E$aR; AReB[rows, , ] <- E$AR; ARcB[rows, , ] <- E$AR }
      else { R0B[rows] <- E0$R; aRcB[rows, ] <- E0$aR; ARcB[rows, , ] <- E0$AR }   # aReB/AReB stay 0 (frozen)
      if (length(lamDir) || .hasCens) {
        s <- .byId[[as.character(.idCode[i])]]; obs <- s[s$EVID == 0, , drop = FALSE]
        if (length(lamDir)) {                          # DV-transform chain (estimated lambda)
          dvSensB[rows, lamDir] <- .foceiAnalyticDvSensLambda(obs$DV, E$trans)
          dvSens2B[rows, lamDir] <- .foceiAnalyticDvSensLambda2(obs$DV, E$trans)
        }
        if (.hasCens) {
          censB[rows] <- if (is.null(obs$CENS)) 0L else as.integer(obs$CENS)
          .lim <- if (is.null(obs$LIMIT)) rep(NA_real_, length(rows)) else as.numeric(obs$LIMIT)
          limB[rows] <- .foceiAnalyticTbsY(.lim, E$trans)
        }
      }
      ehatB[i, ] <- eta0list[[i]]
    }
    nom <- omd$nom
    dOiC <- array(0, c(neta, neta, max(nom, 1L)))
    if (nom > 0L) for (k in seq_len(nom)) dOiC[, , k] <- omd$dOi[[k]]
    d2OiC <- array(0, c(neta, neta, max(nom * nom, 1L)))
    if (nom > 0L) for (aa in seq_len(nom)) for (bb in seq_len(nom)) d2OiC[, , (aa - 1L) * nom + bb] <- omd$d2Oi[[aa]][[bb]]
    d2LD <- if (nom > 0L) omd$d2LD else matrix(0, 1, 1)
    ncores <- tryCatch(as.integer(rxode2::getRxThreads()), error = function(e) 1L)
    if (length(ncores) != 1L || is.na(ncores) || ncores < 1L) ncores <- 1L
    R <- tryCatch(foceiRAllFoceFR_(aB, AB, AthB, aReB, aRcB, AReB, ARcB, dvSensB, dvSens2B,
                                   as.integer(censB), as.numeric(limB), fB, yB, R0B, ehatB, as.integer(off),
                                   Oi, dOiC, d2OiC, d2LD, neta, ndirCov, ndirP, nom, as.integer(dirP), ncores),
                  error = function(e) NULL)
    if (is.null(R) || !all(is.finite(R))) return(NULL)
    return(R)
  }
  # FOCEI: per-subject FD3 (3rd-order Shi) solves collected in R, then ONE OpenMP C++ call
  # (foceiRAllFR_) sums the observed information over subjects -- no per-subject R<->C++ round-trip.
  nsub <- length(ids); Elist <- vector("list", nsub); nobsAll <- integer(nsub)
  # BATCHED (f,R) FOCEI solves: ONE SolveAllFD3 gives f/a/A/Ath AND R/aR/AR/AthR for all subjects
  # (withR=TRUE), instead of the per-subject Shi.  Per-subject Shi is the fallback.
  .obsAll <- lapply(seq_len(nsub), function(i) { .s <- .byId[[as.character(.idCode[i])]]
    if (is.null(.s) || nrow(.s) == 0L) NULL else .s[.s$EVID == 0, , drop = FALSE] })
  if (any(vapply(.obsAll, is.null, logical(1L)))) return(NULL)
  .obsT <- lapply(.obsAll, function(.o) .o$TIME)
  .batch <- !nzchar(Sys.getenv("FOCEI_NO_FD3_BATCH"))
  .EsAll <- if (.batch) .foceiAnalyticSolveAllFD3(am, th, ebes, .idCode, data, .obsT, tol = solveTol, withR = TRUE) else NULL
  if (.batch && is.null(.EsAll)) .batch <- FALSE
  for (i in seq_len(nsub)) {
    s <- .byId[[as.character(.idCode[i])]]; obs <- .obsAll[[i]]
    E <- if (.batch) .EsAll[[i]] else .foceiAnalyticSolveSubjectFD3(am, c(th, setNames(ebes[i, ], etav)), s, obs$TIME, tol = solveTol, withR = TRUE)
    if (is.null(E)) return(NULL)
    if (isTRUE(ef$canVanish)) { .fa <- abs(E$f)
      if (any(!is.finite(.fa)) || min(.fa) < 1e-6 * max(.fa))
        return(.foceiAnalyticFallback("pure proportional error with a near-zero model prediction")) }
    E$y <- .foceiAnalyticTbsY(obs$DV, E$trans)
    Elist[[i]] <- E; nobsAll[i] <- length(E$f)
  }
  totObs <- sum(nobsAll); off <- c(0L, cumsum(nobsAll)); nd2 <- ndirCov * ndirCov
  aB <- matrix(0, totObs, ndirCov); aRB <- matrix(0, totObs, ndirCov)
  AB <- array(0, c(totObs, ndirCov, ndirCov)); ARB <- array(0, c(totObs, ndirCov, ndirCov))
  AthB <- array(0, c(totObs, neta, nd2)); AthRB <- array(0, c(totObs, neta, nd2))
  fB <- numeric(totObs); yB <- numeric(totObs); RB <- numeric(totObs); ehatB <- matrix(0, nsub, neta)
  dvSensB <- if (length(lamDir)) matrix(0, totObs, ndirCov) else matrix(0, totObs, 0L)
  dvSens2B <- dvSensB
  # censored (M2/M3/M4): per-obs CENS + transformed LIMIT drive the censored score partials
  # (the determinant stays Gauss-Newton, matching censOption="gauss"); empty when no censoring.
  .hasCens <- (!is.null(data$CENS) && any(data$CENS != 0, na.rm = TRUE)) ||
    (!is.null(data$LIMIT) && any(is.finite(data$LIMIT)))
  censB <- if (.hasCens) integer(totObs) else integer(0)
  limB <- if (.hasCens) rep(NA_real_, totObs) else numeric(0)
  for (i in seq_len(nsub)) {
    E <- Elist[[i]]; no <- nobsAll[i]; rows <- (off[i] + 1L):off[i + 1L]
    aB[rows, ] <- E$a; aRB[rows, ] <- E$aR; AB[rows, , ] <- E$A; ARB[rows, , ] <- E$AR
    AthB[rows, , ] <- array(E$Ath, c(no, neta, nd2)); AthRB[rows, , ] <- array(E$AthR, c(no, neta, nd2))
    fB[rows] <- E$f; yB[rows] <- E$y; RB[rows] <- E$R; ehatB[i, ] <- ebes[i, ]
    if (length(lamDir) || .hasCens) {
      s <- .byId[[as.character(.idCode[i])]]; obs <- s[s$EVID == 0, , drop = FALSE]
      if (length(lamDir)) {                            # DV-transform chain (estimated lambda)
        dvSensB[rows, lamDir] <- .foceiAnalyticDvSensLambda(obs$DV, E$trans)
        dvSens2B[rows, lamDir] <- .foceiAnalyticDvSensLambda2(obs$DV, E$trans)
      }
      if (.hasCens) {
        censB[rows] <- if (is.null(obs$CENS)) 0L else as.integer(obs$CENS)
        .lim <- if (is.null(obs$LIMIT)) rep(NA_real_, length(rows)) else as.numeric(obs$LIMIT)
        limB[rows] <- .foceiAnalyticTbsY(.lim, E$trans)   # transform the censoring bound like the DV
      }
    }
  }
  nom <- omd$nom
  dOiC <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiC[, , k] <- omd$dOi[[k]]
  d2OiC <- array(0, c(neta, neta, max(nom * nom, 1L)))
  if (nom > 0L) for (aa in seq_len(nom)) for (bb in seq_len(nom)) d2OiC[, , (aa - 1L) * nom + bb] <- omd$d2Oi[[aa]][[bb]]
  d2LD <- if (nom > 0L) omd$d2LD else matrix(0, 1, 1)
  ncores <- tryCatch(as.integer(rxode2::getRxThreads()), error = function(e) 1L)
  if (length(ncores) != 1L || is.na(ncores) || ncores < 1L) ncores <- 1L
  R <- tryCatch(foceiRAllFR_(aB, AB, AthB, aRB, ARB, AthRB, dvSensB, dvSens2B, as.integer(censB), as.numeric(limB),
                             fB, yB, RB, ehatB, as.integer(off),
                             Oi, dOiC, d2OiC, d2LD, neta, ndirCov, ndirP, nom, as.integer(dirP), ncores),
                error = function(e) NULL)
  if (is.null(R) || !all(is.finite(R))) return(NULL)
  R
}

#' Sum the per-subject analytic observed-information R (theta+sigma+Omega) over
#' subjects; `NULL` on any build/solve/non-finite failure.  `startedEnv` (when
#' given) is flagged `.analyticStarted` before the first solve so the C++ hook
#' skips its FD fallback rather than run it on the replaced global solve.
#' @noRd
.foceiAnalyticAssembleR <- function(ui, th, ebes, ids, data, Om, ef, neta, nth, nsg, omd,
                                    dirs, dirTh, ndir, startedEnv = NULL, solveTol = 1e-10,
                                    iovDirScale = NULL, etaScale = NULL, interaction = 1L,
                                    foceType = 0L) {
  # Flag BEFORE building the augmented model: .foceiAnalyticAugModelDirs loads a new
  # compiled model, which calls rxSolveFree() on the fit's global solve.  If the build
  # then declines (is.null(am)/ndir mismatch), the C++ hook must still restore the freed
  # solve before the FD fallback -- otherwise the FD sandwich solves against freed memory.
  if (!is.null(startedEnv)) assign(".analyticStarted", TRUE, startedEnv)
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
  # Batch the 3rd-order solve across ALL subjects (FOCEI *and* FOCE, no IOV rescale) so both take
  # the IDENTICAL method and both get the speedup (1 + 2*neta population solves vs the per-subject
  # Shi's O(nsub*neta)).  CORRECTED FOCE freezes R0 at the eta=0 population solve (batched) and
  # re-solves each EBE (per-subject Newton, censoring-aware) to S_FOCE=0; that eta0 matrix feeds
  # the batched FD3.  foce+ keeps the live R (no eta=0 solve).  Per-subject Shi is the fallback.
  .obsAll <- lapply(seq_along(ids), function(i) { .s <- .byId[[as.character(.idCode[i])]]
    if (is.null(.s) || nrow(.s) == 0L) NULL else .s[.s$EVID == 0, , drop = FALSE] })
  if (any(vapply(.obsAll, is.null, logical(1L)))) return(NULL)   # unmatched subject -> caller FD
  .obsT <- lapply(.obsAll, function(.o) .o$TIME)
  .needF0 <- .foce && identical(as.integer(foceType), 0L) && isTRUE(ef$dependsF0)
  E0List <- vector("list", length(ids))
  if (.needF0) {                                       # eta=0 population solve (frozen R0), batched
    E0List <- .foceiAnalyticSolveAll(am, th, matrix(0, length(ids), neta), .idCode, data, .obsT, solveTol)
    if (is.null(E0List)) return(NULL)
  }
  eta0Mat <- ebes
  if (.foce) for (i in seq_along(ids)) {              # FOCE EBE re-solve (per-subject Newton)
    .o <- .obsAll[[i]]
    .e0 <- .foceiAnalyticFoceEbe(am, th, ebes[i, ], .byId[[as.character(.idCode[i])]], .o$TIME, .o$DV, etav,
                                 if (is.null(E0List[[i]])) NULL else E0List[[i]]$R, Oi, neta, solveTol,
                                 foceType = foceType, cens = .o$CENS, limit = .o$LIMIT)
    if (is.null(.e0)) return(NULL)
    eta0Mat[i, ] <- .e0
  }
  .batch <- !rescale && !nzchar(Sys.getenv("FOCEI_NO_FD3_BATCH"))
  .EsAll <- NULL
  if (.batch) {
    .EsAll <- .foceiAnalyticSolveAllFD3(am, th, eta0Mat, .idCode, data, .obsT, tol = solveTol, withR = FALSE)
    if (is.null(.EsAll)) .batch <- FALSE              # batched solve failed -> per-subject fallback
  }
  for (i in seq_along(ids)) {
    s <- .byId[[as.character(.idCode[i])]]
    obs <- .obsAll[[i]]
    eta0 <- eta0Mat[i, ]
    E0 <- E0List[[i]]
    E <- if (.batch) .EsAll[[i]] else .foceiAnalyticSolveSubjectFD3(am, c(th, setNames(eta0, etav)), s, obs$TIME, tol = solveTol)
    if (is.null(E)) return(NULL)                       # solve failure -> caller falls back to FD
    # near-zero-prediction guard for a vanishing residual variance (pure proportional
    # R = sp^2 f^2 -> 0): the 1/R observed-information terms blow up, so drop to FD.
    if (isTRUE(ef$canVanish)) {
      .fa <- abs(E$f)
      if (any(!is.finite(.fa)) || min(.fa) < 1e-6 * max(.fa))
        return(.foceiAnalyticFallback("pure proportional error with a near-zero model prediction"))
    }
    E$y <- .foceiAnalyticTbsY(obs$DV, E$trans)
    ehat <- eta0
    if (rescale) {
      E$a <- sweep(E$a, 2, iovDirScale, `*`)           # a_B = a_A / w on occasion directions
      for (d1 in seq_len(ndir)) for (d2 in seq_len(ndir))
        E$A[, d1, d2] <- E$A[, d1, d2] * iovDirScale[d1] * iovDirScale[d2]
      for (d1 in seq_len(neta)) for (d2 in seq_len(ndir)) for (d3 in seq_len(ndir))  # Ath is [obs,neta,ndir,ndir]
        E$Ath[, d1, d2, d3] <- E$Ath[, d1, d2, d3] * iovDirScale[d1] * iovDirScale[d2] * iovDirScale[d3]
      if (!is.null(E0)) {                              # population sensitivities share the rescaling
        E0$a <- sweep(E0$a, 2, iovDirScale, `*`)
        for (d1 in seq_len(ndir)) for (d2 in seq_len(ndir))
          E0$A[, d1, d2] <- E0$A[, d1, d2] * iovDirScale[d1] * iovDirScale[d2]
      }
      ehat <- ehat * etaScale                          # xi_hat = w * eta_hat on occasion etas
    }
    Ri <- tryCatch(.foceiAnalyticSubjectR(E, ehat, Om, ef, neta, nth, nsg, ef$sgVar, omd,
                                          ndir = ndir, dirTh = dirTh, Oi = Oi, interaction = interaction,
                                          E0 = E0, foceType = foceType),
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
    if (!.hasRxSens()) return(.foceiAnalyticFallback("an rxode2 without symbolic sensitivities (needs rxExpandSens2_ + symengine)"))
    # bounded thetas are estimated on a transformed scale and the pre-final Jacobian
    # hook corrects env$cov to the natural scale; an analytic install would overwrite
    # that with the internal-scale cov -> bow out to the (Jacobian-correct) FD path.
    if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L)
      return(.foceiAnalyticFallback("a bounded parameter transform"))
    # dose-history functions (tad/podo/tafd/tlast/tfirst/dosenum) are functions of
    # time and the dose record only -- they carry no eta/theta dependence, so rxode2
    # differentiates them to zero (.rxToSEDualVarFunction) and they no longer force
    # the finite-difference fallback.
    # scope: conditional methods only (FOCEI and both FOCE variance modes) with a
    # single Gaussian endpoint; FO/FOI and anything else -> FD.
    if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE))))
      return(.foceiAnalyticFallback("the FO/FOI method"))
    # linCmt() has no symbolic state sensitivities for the augmented model
    if (isTRUE(any(ui$predDf$linCmt)))
      return(.foceiAnalyticFallback("a linCmt() model"))
    interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))                   # 1 FOCEI / 0 FOCE
    # foceType picks the FOCE variance mode (0 "nonmem" frozen R0, 1 "foce+" live R);
    # it only matters when interaction=0 (FOCEI always uses the live conditional R).
    foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(ui, "foceType", 0L)) else 0L
    if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L)
      return(.foceiAnalyticFallback("adaptive Gaussian quadrature (nAGQ > 1)"))
    ef <- .foceiAnalyticErrFull(ui)
    if (is.null(ef)) return(NULL)
    # Estimated boxCox/yeoJohnson lambda: both the FOCEI and FOCE (nonmem / foce+)
    # observed-information cov carry the DV-transform 2nd-order chain (dy'/dlambda residual
    # split in the rho DATA terms + d2y'/dlambda2 in the lambda-lambda block).
    # Both FOCEI and FOCE use the general (f,R) cov path for a general (foceiOnly)
    # variance structure; add/prop keeps the fast symbolic assembly.

    ini <- ui$iniDf
    .map <- .foceiEtaThetaMap(ui)
    etaNames <- .map$etaNames; neta <- length(etaNames)
    if (neta == 0L) return(.foceiAnalyticFallback("a model with no random effects"))
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
          !identical(rxode2::rxGetControl(ui, "iovXform", "sd"), "sd"))
      return(.foceiAnalyticFallback("IOV with a non-SD parameterization"))
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
    ef$ev <- local({ v <- .valc; function(expr, f, y, f0 = f) eval(expr, c(list(f = f, y = y, f0 = f0), as.list(v))) })

    # scope: fixed structural thetas break the eta<->theta indexing
    if (any(.iniIsFixed(ini, thetaForEta)))
      return(.foceiAnalyticFallback("a fixed mu-referenced structural parameter"))
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
      if (length(oe) == 0L || !is.finite(w) || w <= 0)                # unidentified/zero IOV -> FD
        return(.foceiAnalyticFallback("an unidentified or zero IOV variance"))
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
    .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta, iovVars,
                                     sharedEta = unname(.foceiEtaOccurrence(ui) > 1L))
    if (is.null(.dir)) return(NULL)
    thStruct <- .dir$thStruct
    dirs <- .dir$dirs; dirTh <- .dir$dirTh; ndir <- .dir$ndir; nth <- .dir$nth
    # occasion-eta directions carry a[,occ] = w * a_B; rescale by 1/w to the variance-w^2 basis
    iovDirScale <- rep(1, ndir)
    for (g in seq_along(iovVars)) iovDirScale[occEta[[g]]] <- 1 / iovW[g]

    nsg <- length(ef$sgVar)
    # om-param order matches .omegaVarCovDeriv: ordinary Omega elements THEN the IOV
    # shared-variance params (reported on the SD scale as `v`, matching the theta name).
    onm <- etaNames                                            # Omega named by the eta (om.eta.cl)
    # R matrix param order = dirP = [dirTh (structural thetas THEN estimated-lambda) | dirSg
    # (sigmas)] then Omega.  An estimated boxCox/yeoJohnson lambda is a theta-like DIRECTION
    # (already in thStruct/dirTh), so it must NOT be re-listed as a sigma -- use the direction
    # set's sigma names (.dir$sgName, lambda-excluded), not ef$sgName, to avoid a duplicate
    # name that would make fullNm longer than Rfull.
    fullNm <- c(thStruct, .dir$sgName, .foceiOmegaCovNames(pairs, onm), iovVars)  # full natural-scale order

    # op_focei cov-params (non-skipped structural + residual thetas) must all live in
    # the full natural-scale cov (structural thetas + estimated lambda in `thStruct`, residual
    # sigmas in `.dir$sgName`); the R.0 / covFull=FALSE block is taken from it by name below.
    if (!all(covParams %in% fullNm)) return(NULL)

    th <- setNames(as.numeric(thVals[thNames]), paste0("THETA_", seq_along(thNames), "_"))
    etav <- paste0("ETA_", seq_len(neta), "_")
    etaObf <- get("etaObf", e)
    ebes <- as.matrix(etaObf[, paste0("ETA[", seq_len(neta), "]"), drop = FALSE])
    ids  <- etaObf$ID
    data <- get("dataSav", e)
    # Censored data (M2/M3/M4): the analytic observed-information cov carries the censored
    # determinant partials through the general (f,R) path (AssembleRFR), matching the exact
    # censored analytic outer gradient.  Any censored/BLOQ observation (CENS != 0 or a finite
    # LIMIT) routes to that path (and uses the same batched FD3 solve as the uncensored (f,R) cov).
    .hasCensD <- (!is.null(data$CENS) && any(data$CENS != 0, na.rm = TRUE)) ||
      (!is.null(data$LIMIT) && any(is.finite(data$LIMIT)))
    # only the laplace censored determinant stays on the FD cov; gauss (default) is analytic
    if (.hasCensD && as.integer(rxode2::rxGetControl(ui, "censOption", 0L)) == 1L)
      return(.foceiAnalyticFallback("censored observations with censOption='laplace'"))

    # Full natural-scale observed-information R (theta + sigma + Omega), summed over
    # subjects.  startedEnv=e flags `.analyticStarted` before the augmented solve so
    # the C++ hook skips its finite-difference fallback if the assembly then fails.
    # FOCEI with a general (non-add/prop / multi-endpoint) variance uses the (f,R) cov
    # (any structure, sigmas as directions); add/prop FOCEI keeps the fast symbolic
    # assembly until the (f,R) cov is ported to C++ (the R version is correct but slow).
    # FOCE and IOV keep the symbolic add/prop assembly.
    # censored FOCEI must use the general (f,R) path (the fast add/prop assembler has no
    # censored partials); it also carries any general/estimated-lambda variance.
    Rfull <- if (length(iovVars) == 0L && (isTRUE(ef$foceiOnly) || .hasCensD))
      .foceiAnalyticAssembleRFR(ui, th, ebes, ids, data, Om, ef, neta, length(.dir$dirP), .dir$dirP, omd,
                                dirsCov = .dir$dirsCov, ndirCov = .dir$ndirCov,
                                startedEnv = e, solveTol = .foceiAnalyticSolveTol(ui),
                                interaction = interaction, foceType = foceType, lamDir = .dir$lamDir)
    else .foceiAnalyticAssembleR(ui, th, ebes, ids, data, Om, ef, neta, nth, nsg, omd,
                                 dirs = dirs, dirTh = dirTh, ndir = ndir,
                                 startedEnv = e, solveTol = .foceiAnalyticSolveTol(ui),
                                 iovDirScale = iovDirScale, etaScale = etaScale,
                                 interaction = interaction, foceType = foceType)
    if (is.null(Rfull)) return(NULL)
    dimnames(Rfull) <- list(fullNm, fullNm)

    # full natural-scale cov -> fit$cov; theta SEs flow through the native path
    .covNat <- tryCatch(solve(Rfull), error = function(e) NULL)
    if (is.null(.covNat)) return(NULL)
    dimnames(.covNat) <- list(fullNm, fullNm)
    # the op_focei cov-param block (structural + residual thetas, by name) decides
    # success; the native covR = Rinv path (issue #666 fix) reports these SEs directly.
    .ct <- .covNat[covParams, covParams, drop = FALSE]
    .R0 <- tryCatch(solve(.ct), error = function(e) NULL)
    if (is.null(.R0)) return(NULL)   # inversion failed -> FD, and do NOT leave a stale .analyticCov
    assign(".analyticCov", .covNat, envir = e)           # stash only after the deciding inversion
    # covFull=FALSE installs the non-skipped theta block (structural + residual),
    # matching the covType="fd" shape (skipCov drops only fixed/IOV/mixProb thetas).
    assign(".analyticThetaNames", covParams, envir = e)
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
  # Distribution guard (Gaussian-only).  The analytic gradient AND covariance assemble a
  # conditional-GAUSSIAN observed information (rho = 0.5*((y-f)^2/R + log R)).  A non-normal
  # endpoint (t, cauchy, poisson, binomial, ordinal, ...) is built on the SAME norm-form
  # rx_pred_/rx_r_ as a Gaussian fit -- its real (heavier-tailed / discrete) density lives in a
  # SEPARATE ..Llik model the augmented sensitivity model never touches.  Such an endpoint that
  # still carries an add()/prop() location scale would slip past the err-based checks below and
  # be assembled as Gaussian (silently wrong).  The gradient is already shielded upstream by the
  # needOptimHess -> fast=FALSE downgrade (focei.R), but covType="analytic" has no such
  # downgrade, so gate it here.  "norm" and its explicit-Gaussian spelling "dnorm" (used for the
  # Gaussian endpoints of a mixed model) stay in scope; everything else falls back to FD.
  .dist <- tryCatch(as.character(ui$predDfFocei$distribution), error = function(e) NULL)
  if (!length(.dist)) .dist <- tryCatch(as.character(ui$predDf$distribution), error = function(e) character(0))
  if (length(.dist) && !all(.dist %in% c("norm", "dnorm")))
    return(.foceiAnalyticFallback("a non-normal likelihood endpoint (t/cauchy/count/ordinal); the analytic path is Gaussian-only"))
  # Multiple modeled endpoints: rx_pred_ and rx_r_ are single dvid-conditional
  # expressions that already select the right endpoint per observation when solved
  # against the dataset, so the (f,R) path handles them -- but the single-endpoint
  # symbolic add/prop machinery (one rx_pred_) does not, so force the general path.
  .multiEndpoint <- !is.null(ui$predDf) && nrow(ui$predDf) > 1L
  ini <- ui$iniDf
  # Both-sides transforms: rx_pred_ is already the transformed prediction and rx_r_ the
  # transformed-scale variance, so the analytic path only transforms the DV to that scale
  # (y' = tbs(DV), via the C++ _powerD in the assemblers).  Exact for a STATIC transform
  # (lnorm, or a FIXED boxCox/yeoJohnson lambda): the transform Jacobian -2 log|dy'/dy| is
  # then constant (or linear in a fixed lambda) and drops out of the observed information.
  # An ESTIMATED lambda additionally makes the DV move with the parameter (a dy'/dlambda
  # chain in the residual): the GRADIENT carries it (lambda is a theta-like direction with
  # df'/dlambda plus the residual DV chain and the -2 log|J| Jacobian term); the cov path
  # gates estimated lambda separately (its 2nd-order lambda terms are not yet ported).
  .trans <- as.character(ui$predDf$transform)
  .estLam <- FALSE
  if (!all(.trans == "untransformed")) {
    .lamRows <- ini[!is.na(ini$err) & ini$err %in% c("boxCox", "yeoJohnson"), , drop = FALSE]
    .estLam <- any(!.lamRows$fix)
    # Multiple estimated per-endpoint lambdas need an endpoint->lambda DV mapping that is not
    # wired (the gradient gates this at foceiGradAnalytic.R; without the same gate the covariance
    # would recycle one endpoint's dy'/dlambda chain across every lambda column -- silently
    # wrong).  A single estimated lambda is the ported case; more than one falls back to FD.
    if (sum(!.lamRows$fix) > 1L)
      return(.foceiAnalyticFallback("multiple estimated boxCox/yeoJohnson lambdas (endpoint->lambda mapping not yet wired)"))
  }
  er <- ini[!is.na(ini$err), , drop = FALSE]
  if (nrow(er) == 0L) {
    # All residual-error parameters are fixed: the fit's working iniDf drops
    # them (they are baked into the compiled rx_r_ as constants), yet the model
    # still HAS a residual (predDf carries its errType).  Rather than fall back
    # to FD -- which re-solves the ODE per outer parameter (all omegas here,
    # which do not even enter the ODE) -- run the general (f,R) path with NO
    # free sigma direction (nsg = 0): rx_r_ is read from the solve, the omega
    # gradient uses the precomputed Omega derivatives, and R's constancy needs
    # no re-solve.  Only the plain conditional-Gaussian error types are ported;
    # anything else stays on FD.
    .et <- as.character(ui$predDf$errType)
    .tok <- unique(trimws(unlist(strsplit(.et, "[^A-Za-z0-9]+"))))
    .tok <- .tok[nzchar(.tok)]
    .known <- c("add", "prop", "pow", "combined1", "combined2")
    if (length(.tok) == 0L || !all(.tok %in% .known) ||
          !all(.trans == "untransformed"))
      return(.foceiAnalyticFallback("a model with no estimated residual error"))
    .hasAddFloor <- any(.tok %in% c("add", "combined1", "combined2"))
    return(list(sgVar = character(0), sgName = character(0), sc = NULL,
                per = NULL, pair = NULL, foce = NULL, focePlus = NULL,
                foceiOnly = TRUE, canVanish = !.hasAddFloor,
                dependsF0 = !all(.tok == "add"), estLam = .estLam,
                ev = function(e, f, y, f0 = f) NULL))
  }
  sgNameAll <- er$name                               # ALL error params (excluded from directions)
  # pure proportional / power error (no additive floor) vanishes as f -> 0, making the
  # 1/R observed-information terms blow up near zero predictions; the assembly guards it.
  # lnorm/logitNorm/probitNorm are additive on their transformed (log / logit / probit) scale
  # (R = sd^2 constant), so the variance never vanishes even though the transformed prediction
  # crosses zero (e.g. logit(cp) = 0 at the bound midpoint) -- the canVanish guard, meant for
  # pure proportional/power error whose R -> 0 with f, must not fire for them.
  canVanish <- !any(er$err %in% c("add", "lnorm", "logitNorm", "probitNorm"))
  # model-declared addProp wins; the control applies only when the model says "default"
  addPr <- as.character(ui$predDf$addProp)
  if (length(addPr) != 1L || is.na(addPr) || addPr == "default") {
    addPr <- tryCatch(rxode2::rxGetControl(ui, "addProp", "combined2"), error = function(e) "combined2")
  }
  # The general (f,R) path reads R and dR/dsigma from the model's own rx_r_ via the solve, so ANY
  # conditional-Gaussian variance structure is in scope -- including combined1 (sa+sp*f)^2, which
  # rx_r_ carries directly ((add + pred*prop)^2, the weight squared).  The symbolic add/prop path
  # below is only the FAST route for the plain untransformed add/prop/combined2 case; everything
  # else returns a minimal `ef` (foceiOnly=TRUE) and runs the general (f,R) assembler
  # (.foceiAnalyticAssembleRFR), which handles BOTH FOCEI and FOCE (interaction 0/1).
  # a both-sides transform is never the plain-scale symbolic add/prop machinery: force the
  # general (f,R) path (rx_pred_/rx_r_ carry the transform; only the DV is retransformed).
  .isAddProp <- !.multiEndpoint && all(.trans == "untransformed") &&
    all(er$err %in% c("add", "prop")) && !identical(addPr, "combined1")
  # R0 (FOCE nonmem frozen variance) needs the eta=0 population solve only when R depends
  # on the prediction (any non-additive error term); pure additive R is constant.
  .dependsF0 <- !all(er$err == "add")
  if (!.isAddProp)
    return(list(sgVar = character(0), sgName = sgNameAll, sc = NULL, per = NULL, pair = NULL,
                foce = NULL, focePlus = NULL, foceiOnly = TRUE, canVanish = canVanish,
                dependsF0 = .dependsF0, estLam = .estLam, ev = function(e, f, y, f0 = f) NULL))
  addN <- er$name[er$err == "add"]; propN <- er$name[er$err == "prop"]
  hasA <- length(addN) == 1L; hasP <- length(propN) == 1L
  if (!hasA && !hasP)
    return(.foceiAnalyticFallback("a model with no additive or proportional residual error"))
  Rstr <- if (hasA && hasP) "sa^2+sp^2*f^2" else if (hasP) "sp^2*f^2" else "sa^2"
  Rq <- parse(text = Rstr)[[1]]
  rhoE <- bquote(0.5 * ((y - f)^2 / .(Rq) + log(.(Rq))))
  pE <- bquote(1 / .(Rq) + 0.5 * (.(stats::D(Rq, "f")) / .(Rq))^2)
  # FOCE (interaction=0) pieces: the inner problem drops dR/deta, so its gradient
  # coefficient is the least-squares part only, qE = drho/df = -(y-f)/R
  # (q'=dq/df, q''=d2q/df2 build the FOCE inner Hessian Hf and its 3-tensor); the FOCE
  # Laplace determinant curvature is pFE = 1/R (no 0.5*(R'/R)^2 interaction term).
  # R here is written in the live prediction f, so D(, "f") carries the dR/df chain:
  # these are the "foce+" (live conditional R) pieces; the "nonmem" frozen-R0 variants
  # are built below with the separate f0 symbol.
  qE <- bquote(-(y - f) / .(Rq))
  pFE <- bquote(1 / .(Rq))
  DD <- function(e, ...) { for (v in c(...)) e <- stats::D(e, v); e }
  sgVar <- c(if (hasA) "sa", if (hasP) "sp"); sgName <- c(if (hasA) addN, if (hasP) propN)
  val <- setNames(c(if (hasA) er$est[er$name == addN], if (hasP) er$est[er$name == propN]), sgVar)
  sc <- list(r1 = DD(rhoE, "f"), r2 = DD(rhoE, "f", "f"), r3 = DD(rhoE, "f", "f", "f"),
             p = pE, p1 = DD(pE, "f"), p2 = DD(pE, "f", "f"),
             q0 = qE, q1 = DD(qE, "f"), q2 = DD(qE, "f", "f"),
             pF = pFE, pF1 = DD(pFE, "f"), pF2 = DD(pFE, "f", "f"))
  per <- list(); for (s in sgVar) per[[s]] <- list(rf = DD(rhoE, "f", s), rff = DD(rhoE, "f", "f", s),
    ps = DD(pE, s), pf = DD(pE, "f", s),
    qs = DD(qE, s), qsf = DD(qE, "f", s), psF = DD(pFE, s), pfF = DD(pFE, "f", s),
    rs = DD(rhoE, s))   # pure d(rho)/d(sigma): first-derivative (outer-gradient) sigma term
  pair <- list(); for (i in seq_along(sgVar)) for (j in i:length(sgVar)) { a <- sgVar[i]; b <- sgVar[j]
    pair[[paste0(a, b)]] <- list(rss = DD(rhoE, a, b), rfss = DD(rhoE, "f", a, b), pss = DD(pE, a, b),
      qss = DD(qE, a, b), pssF = DD(pFE, a, b)) }
  # CORRECTED FOCE variance R0 = R at the eta=0 POPULATION prediction (symbol `f0`),
  # while the residual/numerator keeps the eta-hat prediction `f` (the estimator freezes
  # R at eta=0: getPopR/likInner0).  Because f0 is a separate symbol, D(.,"f") treats R0
  # as constant (so q1=1/R0, q2=0, pF1=0, r3=0), and D(.,f0) gives the theta-chain (R0
  # varies with theta through the population prediction) supplied via a0/A0 in the
  # subject assembly.  Additive R0=sa^2 has no f0 -> all f0 fields vanish (== FOCEI-add).
  R0q <- parse(text = gsub("\\bf\\b", "f0", Rstr))[[1]]
  rho0E <- bquote(0.5 * ((y - f)^2 / .(R0q) + log(.(R0q))))
  q0E <- bquote(-(y - f) / .(R0q)); pF0E <- bquote(1 / .(R0q))
  scF <- list(r1 = DD(rho0E, "f"), r2 = DD(rho0E, "f", "f"), r3 = DD(rho0E, "f", "f", "f"),
              q0 = q0E, q1 = DD(q0E, "f"), q2 = DD(q0E, "f", "f"),
              pF = pF0E, pF1 = DD(pF0E, "f"), pF2 = DD(pF0E, "f", "f"))
  perF <- list(); for (s in sgVar) perF[[s]] <- list(rf = DD(rho0E, "f", s), rff = DD(rho0E, "f", "f", s),
    qs = DD(q0E, s), qsf = DD(q0E, "f", s), psF = DD(pF0E, s), pfF = DD(pF0E, "f", s),
    rs = DD(rho0E, s))   # pure d(rho)/d(sigma) for the FOCE first-derivative sigma term
  pairF <- list(); for (i in seq_along(sgVar)) for (j in i:length(sgVar)) { a <- sgVar[i]; b <- sgVar[j]
    pairF[[paste0(a, b)]] <- list(rss = DD(rho0E, a, b), rfss = DD(rho0E, "f", a, b),
      qss = DD(q0E, a, b), pssF = DD(pF0E, a, b)) }
  f0F <- list(qf0 = DD(q0E, "f0"), qff0 = DD(q0E, "f", "f0"), qf0f0 = DD(q0E, "f0", "f0"),
              pFf0 = DD(pF0E, "f0"), pFf0f0 = DD(pF0E, "f0", "f0"),
              rhof0 = DD(rho0E, "f0"), rhof0f0 = DD(rho0E, "f0", "f0"))
  perf0F <- list(); for (s in sgVar) perf0F[[s]] <- list(qf0s = DD(q0E, "f0", s),
    pFf0s = DD(pF0E, "f0", s), rhof0s = DD(rho0E, "f0", s))
  foce <- list(sc = scF, per = perF, pair = pairF, f0 = f0F, perf0 = perf0F, dependsF0 = hasP)
  # "foce+" (foceType=1, live conditional R): the same truncated inner gradient and
  # determinant, but R follows the eta-hat prediction, so the live-R sc/per/pair pieces
  # apply as-is and there is no f0 population solve/chain.
  focePlus <- list(sc = sc, per = per, pair = pair, dependsF0 = FALSE)
  list(sgVar = sgVar, sgName = sgName, sc = sc, per = per, pair = pair, foce = foce,
       focePlus = focePlus, foceiOnly = FALSE, canVanish = canVanish, dependsF0 = .dependsF0,
       estLam = .estLam,
       ev = function(e, f, y, f0 = f) eval(e, c(list(f = f, y = y, f0 = f0), as.list(val))))
}

#' Augmented rxode2 model with state + 1st/2nd-order sensitivities and the
#' prediction chain f1/f2 over an arbitrary direction set `dirs` (each ETA_i_ or
#' THETA_j_).  State sensitivities use rxode2's `.rxSens`; the higher-order
#' prediction chain f1/f2 is built here.  O(ndir^2), so it compiles far more
#' cheaply than a 3rd-order model.
#' @return list(augMod, dirs, ndir, st, P2) or `NULL` on failure
#' @noRd
.foceiAnalyticAugCache <- new.env(parent = emptyenv())     # session cache: (model digest | dirs) -> aug model
#' Solve-output column names + index maps for the augmented model (from dirs, P2, sigTh).
#' Precomputed once on `am$cols` at build time so the matrix-extraction readers avoid per-call
#' paste0; recomputed on the fly for a reconstructed `am` (from foceiModel$outer) that predates it.
#' @noRd
.foceiAnalyticCols <- function(dirs, fDirs, P2, P2r, sigTh) {
  sigP2 <- if (length(sigTh)) do.call(rbind, lapply(seq_along(sigTh), function(.a)
    data.frame(a = .a, b = seq_along(sigTh)[seq_along(sigTh) >= .a]))) else NULL
  # f-part (prediction a/A) spans f-directions only (sigma slots are zero-filled by the reader);
  # R-part (aR/AR) spans every direction.  iiF/jjF index the f2 pairs into the FULL dirs vector.
  list(f1 = paste0("rx_f1_", fDirs), fDirIdx = match(fDirs, dirs),
       f2 = paste0("rx_f2_", P2$i, "_", P2$j), iiF = match(P2$i, dirs), jjF = match(P2$j, dirs),
       rvar1 = paste0("rx_rvar1_", dirs), rvar2 = paste0("rx_rvar2_", P2r$i, "_", P2r$j),
       ii = match(P2r$i, dirs), jj = match(P2r$j, dirs),
       rsig = if (length(sigTh)) paste0("rx_rsig_", sigTh, "_") else character(0),
       rsig1 = lapply(sigTh, function(.n) paste0("rx_rsig1_", .n, "_", dirs)),
       rsig2 = if (is.null(sigP2)) character(0) else paste0("rx_rsig2_", sigTh[sigP2$a], "_", sigTh[sigP2$b]),
       sigP2 = sigP2)
}
.foceiAnalyticEtCache <- new.env(parent = emptyenv())      # per-fit translated event table (etTrans reuse)
#' The augmented model's events are IDENTICAL across every solve of a fit (only the theta/eta
#' params vary), so translate them once with etTrans and reuse -- ~40% off each population solve
#' (validated: identical predictions).  Keyed by the model key + a content hash of the data
#' (`digest`, already a dependency), so distinct datasets that happen to share summary statistics
#' cannot collide onto the same translated event table; falls back to raw data.  The hash is far
#' cheaper than the etTrans it guards, so the per-solve lookup stays a net win.
#' @noRd
.foceiAnalyticEvents <- function(am, data) {
  if (is.null(am$key)) return(data)
  .fp <- tryCatch(digest::digest(data), error = function(e) NULL)
  if (is.null(.fp)) return(data)                          # cannot fingerprint -> translate fresh (no cache)
  .ek <- paste0(am$key, "|et|", nrow(data), "|", .fp)
  .et <- get0(.ek, envir = .foceiAnalyticEtCache, inherits = FALSE)
  if (is.null(.et)) {
    .et <- tryCatch(rxode2::etTrans(data, am$augMod), error = function(e) NULL)
    if (is.null(.et)) return(data)
    if (length(ls(.foceiAnalyticEtCache, all.names = TRUE)) >= 256L)    # bound memory in long sessions
      rm(list = ls(.foceiAnalyticEtCache, all.names = TRUE), envir = .foceiAnalyticEtCache)
    assign(.ek, .et, envir = .foceiAnalyticEtCache)
  }
  .et
}
#' Covariate-coefficient directions that can be reconstructed by eta-scaling.
#'
#' A mu-referenced covariate coefficient `b` enters linearly through an eta's parameter
#' (`cl = exp(tcl + eta.cl + b*cov)`), so for a subject-constant covariate every sensitivity
#' of `b` equals the linked eta's scaled by the covariate value: df/db = cov*df/deta,
#' d2f/(db dx) = cov*d2f/(deta dx), d2f/db2 = cov^2*d2f/deta2 (verified exact).  We can therefore
#' skip integrating `b`'s state-sensitivity ODEs and emit its f1/f2/rvar/rsig columns as scaled
#' copies of the linked eta's columns.  Returns a named list `covDir -> list(etaDir, scale)`
#' (scale = the covariate expression in rxode2 syntax) for such thetas present in `fDirs`.
#'
#' Detection reuses rxode2's OWN mu-reference classification rather than re-parsing the model --
#' the same machinery the mu-referenced (muModel="lin"/"irls") family consumes (see #711/#712):
#'   * ui$muRefDataFrame                  : the theta<->eta link.  A well-formed mu-reference has
#'       the eta as a BARE `+eta` in the parameter; a scaled eta (`0.5*eta`) is simply absent, so
#'       such coefficients are auto-excluded with no extra check.
#'   * ui$muRefCovariateDataFrame         : BARE data-column covariate coefficients.
#'   * ui$mu2RefCovariateReplaceDataFrame : ALGEBRAIC covariate-expression coefficients (allometric
#'       `log(WT/70)`, centered `(WT-70)`, ...) -- exactly the set `muRefCovAlg` rewrites; its
#'       `covariate` column is the (equivalent) covariate expression used as the scale.
#' Both covariate frames carry (theta, covariateParameter, covariate); the theta is mapped to its
#' eta via muRefDataFrame.  The ONE structural condition rxode2's frames do not encode is an eta
#' REUSED in another parameter (then df/deta gains a path df/db lacks), so an eta-occurrence guard
#' is applied.  The scale being covariate-and-constant only is enforced downstream by the caller's
#' subject-constant gate.  Empty when the info is unavailable.
#' @noRd
.foceiAnalyticCovDirMap <- function(ui, fDirs) {
  .ini <- ui$iniDf
  # eta name -> ETA_k_ (k = order among diagonal etas, matching .foceiAnalyticDirections' etaDirs)
  .etaRows <- .ini[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2, , drop = FALSE]
  .etaRows <- .etaRows[order(.etaRows$neta1), , drop = FALSE]
  .eta2Dir <- stats::setNames(paste0("ETA_", seq_len(nrow(.etaRows)), "_"), .etaRows$name)
  # id-level theta -> eta (occasion/IOV etas excluded -- IOV is out of analytic scope anyway).
  # May be empty: a covariate can also scale an ETA-LESS structural parameter (handled below), so
  # do NOT bail on an empty eta map -- only bail when there are no covariate coefficients at all.
  .mrd <- tryCatch(ui$muRefDataFrame, error = function(e) NULL)
  .theta2eta <- if (!is.null(.mrd) && is.data.frame(.mrd) && nrow(.mrd) > 0L) {
    .idRef <- .mrd[is.na(.mrd$level) | .mrd$level == "id", , drop = FALSE]
    stats::setNames(as.character(.idRef$eta), as.character(.idRef$theta))
  } else character(0)
  # (theta, coefficient, scale) rows from BOTH rxode2 covariate frames (bare + algebraic)
  .grab <- function(df, scaleCol) {
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0L ||
        !all(c("theta", "covariateParameter", scaleCol) %in% names(df))) return(NULL)
    data.frame(theta = as.character(df$theta), coef = as.character(df$covariateParameter),
               scale = as.character(df[[scaleCol]]), stringsAsFactors = FALSE)
  }
  .cov <- do.call(rbind, list(
    .grab(tryCatch(ui$muRefCovariateDataFrame, error = function(e) NULL), "covariate"),
    .grab(tryCatch(ui$mu2RefCovariateReplaceDataFrame, error = function(e) NULL), "covariate")))
  if (is.null(.cov) || nrow(.cov) == 0L) return(list())
  # A covariate coefficient's sensitivities are those of the parameter position it scales, times the
  # covariate value.  That position is either a random effect (`cl = exp(tcl + eta.cl + b*cov)` ->
  # df/db = cov*df/deta) or, for a covariate on an ETA-LESS parameter (`v = exp(tv + b*cov)`, no
  # eta.v), the structural theta itself (df/db = cov*df/dtheta, since theta and b enter the mu
  # identically -- tv's own direction is already integrated, so b reuses it for free).  Guard the
  # ONE structural condition the rxode2 frames don't encode: the reused position must occur exactly
  # once; a shared eta/theta gains a state path df/db lacks (counts are conservative -> `2L` when
  # lstExpr is unavailable, so uncertain cases fall back to the full symbolic build).
  .eo <- .foceiEtaOccurrence(ui)
  .to <- .foceiNameOccurrence(ui, unique(.cov$theta))
  .res <- list(); .dup <- character(0)
  for (.r in seq_len(nrow(.cov))) {
    .thN <- .cov$theta[.r]
    .etN <- if (.thN %in% names(.theta2eta)) .theta2eta[[.thN]] else NA_character_
    if (!is.na(.etN)) {
      .ed <- .eta2Dir[[.etN]]; if (is.null(.ed) || !(.ed %in% fDirs)) next
      if (!identical(.eo[[.etN]], 1L)) next                      # eta shared/reused -> scaling not exact
    } else {
      # eta-less parameter: reuse the structural theta's own (already integrated) direction
      .snt <- .ini$ntheta[.ini$name == .thN]; if (length(.snt) != 1L || is.na(.snt)) next
      .ed <- paste0("THETA_", .snt, "_"); if (!(.ed %in% fDirs)) next
      if (!identical(.to[[.thN]], 1L)) next                      # structural theta shared -> scaling not exact
    }
    .nt <- .ini$ntheta[.ini$name == .cov$coef[.r]]; if (length(.nt) != 1L || is.na(.nt)) next
    .cd <- paste0("THETA_", .nt, "_"); if (!(.cd %in% fDirs) || .cd %in% .dup) next
    # `etaDir` names the source direction to scale (an eta direction, or a structural-theta
    # direction for an eta-less parameter); the downstream emitter treats both identically.
    .new <- list(etaDir = .ed, scale = .cov$scale[.r])
    if (!is.null(.res[[.cd]])) {
      # same coefficient re-listed identically (bare frame + algebraic frame overlap) -> no-op; a
      # different (etaDir, scale) is a genuine ambiguity (b*cov1 + b*cov2) -> drop, fall back.
      if (identical(.res[[.cd]], .new)) next
      .res[[.cd]] <- NULL; .dup <- c(.dup, .cd); next
    }
    .res[[.cd]] <- .new
  }
  .res
}
.foceiAnalyticAugModelDirs <- function(ui, dirs) {
  .key <- tryCatch(paste0(rxUiGet.foceiModelDigest(list(ui)), "|", paste(dirs, collapse = ","),
                          "|sk", Sys.getenv("FOCEI_NO_SIGMA_SKIP"),     # skip flag -> distinct cached model
                          # subject-constant covariate set -> which covariate directions get eta-scaling reuse
                          "|cc", paste(sort(tryCatch(rxode2::rxGetControl(ui, "foceiConstCovs", NULL),
                                                     error = function(e) NULL)), collapse = ",")),
                   error = function(e) NULL)
  if (!is.null(.key)) { .hit <- get0(.key, envir = .foceiAnalyticAugCache, inherits = FALSE)
    if (!is.null(.hit)) return(.hit) }
  .res <- tryCatch({
    .s <- ui$loadPruneSens
    .st <- rxode2::rxStateOde(.s)
    # matExp()/indLin(): the base ODEs are materialized into the pruned env (via
    # .rxJacobian -> .rxInjectMatExpOdes) with correct RHS, but rxStateOde() can return
    # them reversed (an indLin() state parses as compartment 1).  The augmented model
    # numbers compartments by first d/dt() appearance, so emitting .baseOde in that
    # reversed order misplaces default dosing (the dose lands in the wrong compartment)
    # and the eta sensitivities collapse.  Reorder source-first, matching the inner model
    # (.rxInjectMatExpDdt / rxUiGet.foceiCmtPreModel).
    .mvS <- rxode2::rxModelVars(.s)
    if (is.list(.mvS$indLin) && length(.mvS$indLin) == 4L) {
      .st <- .rxMatExpStateOrder(.st, ls(envir = .s, all.names = TRUE))
    }
    # Sigma-skip (default ON; set FOCEI_NO_SIGMA_SKIP=1 to opt out): the residual-error (sigma) directions have
    # df/dsigma=0, so their f-sensitivity columns (f1/f2) and state-sensitivities are all
    # ZERO -- pure wasted compile (each sigma direction costs a full structural parameter's
    # worth of gcc + states).  Build the prediction chain f1/f2 + state sensitivities over
    # the model-affecting directions only (.fDirs); the variance chain still spans EVERY
    # direction (sigma R-derivatives are algebraic, no state chain).  The reader zero-fills
    # the sigma a/A slots (the assembly already treats a=A=Ath=0 for a sigma direction).
    .erN <- ui$iniDf$ntheta[!is.na(ui$iniDf$err) & !(ui$iniDf$err %in% c("boxCox", "yeoJohnson"))]
    .sigDirs <- if (!nzchar(Sys.getenv("FOCEI_NO_SIGMA_SKIP"))) intersect(dirs, paste0("THETA_", .erN, "_")) else character(0)
    .fDirs <- dirs[!(dirs %in% .sigDirs)]
    # Covariate-coefficient reuse (always on; the subject-constant gate below keeps it exact): a
    # mu-ref covariate direction's sensitivities are the linked eta's scaled by the covariate value,
    # so we build state sensitivities over the model directions .mfDirs only (no covariate-direction
    # ODEs) and emit the covariate columns algebraically below.  The covariate directions stay in the
    # full direction set (.fDirs/dirs/.P2/.P2r/am$cols) so the assembly and FD3 are unchanged -- they
    # read the algebraically-emitted covariate columns like any other.
    .covMap <- .foceiAnalyticCovDirMap(ui, .fDirs)
    # Subject-constant gate: the eta-scaling identity holds only for covariates that are constant
    # within each subject (a time-varying covariate makes d(state)/db != cov*d(state)/deta).  The
    # constant-covariate set is computed from the data and stashed on the control by the fit setup;
    # reuse only a covariate whose scale expression uses exclusively constant covariates (and drop
    # all reuse when the set is unavailable, e.g. standalone use -- the full build is always correct).
    if (length(.covMap)) {
      .constCovs <- tryCatch(rxode2::rxGetControl(ui, "foceiConstCovs", NULL), error = function(e) NULL)
      .covMap <- .covMap[vapply(.covMap, function(.m) {
        .v <- all.vars(parse(text = .m$scale)[[1]]); length(.v) > 0L && all(.v %in% .constCovs)
      }, logical(1))]
    }
    .covDirs <- names(.covMap)
    .mfDirs <- .fDirs[!(.fDirs %in% .covDirs)]      # directions that get integrated state sensitivities
    rxode2::.rxJacobian(.s, c(.st, .mfDirs))
    # 1st-order sensitivities are already expanded for the gradient (free); if
    # unavailable the model is not differentiable, so bail before the 2nd-order build.
    .s1 <- rxode2::.rxSens(.s, .mfDirs)
    if (length(.s1) == 0L) return(NULL)
    .s2 <- rxode2::.rxSens(.s, .mfDirs, .mfDirs)    # 2nd order (model f-directions only): the expensive expansion
    .pred <- get("rx_pred_", .s)
    # residual variance rx_r_ (any structure) with its own 1st/2nd sensitivity chains,
    # exactly like the prediction -- so R and dR/ddir, d2R/ddir2 come from the SOLVE
    # (matching the inner ODE model's d(R)/d(eta)); the assembly then treats the
    # transformed prediction f and the variance R as independent solved quantities.
    .rvar <- tryCatch(get("rx_r_", .s), error = function(e) NULL)
    .Dn <- function(.e, .v) symengine::D(.e, symengine::S(.v))
    .sn1 <- function(.j, ...) symengine::S(paste0("rx__sens_", .j, "_BY_", paste(c(...), collapse = "_BY_"), "__"))
    .toRx <- function(.l) rxode2::rxFromSE(.l)
    # state-sensitivity chains only for model directions .mfDirs (a sigma or covariate-reuse
    # direction contributes no integrated state chain -- sigma has d(state)/dsigma=0, and a
    # covariate direction's columns are emitted by scaling below, not via .g1/.g2).
    .g1 <- function(.ex, .p) { .e <- .Dn(.ex, .p)
      if (.p %in% .mfDirs) for (.j in .st) .e <- .e + .Dn(.ex, .j) * .sn1(.j, .p); .e }
    .g2 <- function(.ex, .p, .q) { .gq <- .g1(.ex, .q); .e <- .Dn(.gq, .p)
      if (.p %in% .mfDirs) for (.k in .st) .e <- .e + .Dn(.gq, .k) * .sn1(.k, .p)
      if (.p %in% .mfDirs && .q %in% .mfDirs) for (.j in .st) .e <- .e + .Dn(.ex, .j) * .sn1(.j, .p, .q); .e }
    # Covariate-direction reuse: map a direction to its model substitute (the linked eta for a
    # covariate direction) + scale; emit a 1st/2nd-order line symbolically when it involves only
    # model directions, else algebraically from the (eta-substituted, reordered) model column
    # scaled by the covariate value(s).  `.emitScaled1` does the 1st-order (single-direction)
    # columns; `.emitPair` the 2nd-order (pair) columns.
    .covSub <- function(.d) { .m <- .covMap[[.d]]; if (is.null(.m)) list(d = .d, s = NULL) else list(d = .m$etaDir, s = .m$scale) }
    .emitScaled1 <- function(.pre, .suf, .cd) paste0(.pre, .cd, .suf, "=(", .covMap[[.cd]]$scale, ")*", .pre, .covMap[[.cd]]$etaDir, .suf)
    .emitPair <- function(.pre, .i, .j, .symFn, .ord) {
      .si <- .covSub(.i); .sj <- .covSub(.j)
      if (is.null(.si$s) && is.null(.sj$s)) return(paste0(.pre, .i, "_", .j, "=", .symFn(.i, .j)))
      .a <- .si$d; .b <- .sj$d
      if (match(.a, .ord) > match(.b, .ord)) { .t <- .a; .a <- .b; .b <- .t }
      paste0(.pre, .i, "_", .j, "=", paste(sprintf("(%s)", c(.si$s, .sj$s)), collapse = "*"), "*", .pre, .a, "_", .b)
    }
    # f2_i_j == f2_j_i, so only the i<=j triangle is differentiated/compiled (~2x off
    # the dominant model-build cost); the reader mirrors A[,i,j]=A[,j,i].
    .P2 <- expand.grid(i = .fDirs, j = .fDirs, stringsAsFactors = FALSE)     # prediction f2: f-directions only
    .P2 <- .P2[match(.P2$i, .fDirs) <= match(.P2$j, .fDirs), , drop = FALSE]
    .P2r <- expand.grid(i = dirs, j = dirs, stringsAsFactors = FALSE)        # variance rvar2: every direction
    .P2r <- .P2r[match(.P2r$i, dirs) <= match(.P2r$j, dirs), , drop = FALSE]
    # model directions first (they define the columns the covariate columns scale), covariate
    # directions after -- rxode2 evaluates assignments in order.
    .fL1 <- c(vapply(.mfDirs, function(.p) paste0("rx_f1_", .p, "=", .toRx(.g1(.pred, .p))), character(1)),
              vapply(.covDirs, function(.cd) .emitScaled1("rx_f1_", "", .cd), character(1)))
    .p2mod <- !(.P2$i %in% .covDirs) & !(.P2$j %in% .covDirs)
    .fL2 <- c(vapply(which(.p2mod), function(.r) .emitPair("rx_f2_", .P2$i[.r], .P2$j[.r], function(.a, .b) .toRx(.g2(.pred, .a, .b)), .fDirs), character(1)),
              vapply(which(!.p2mod), function(.r) .emitPair("rx_f2_", .P2$i[.r], .P2$j[.r], NULL, .fDirs), character(1)))
    # Variance rx_r_ and its 1st/2nd direction sensitivities.  The `rx_..._` naming
    # matters: rxode2 keeps a constant assignment to an `rx_<name>_` variable as a
    # real lhs output column (like the inner model's rx__sens_*), whereas a plainly
    # named `x = <constant>` would be parsed as an initial value and drop out of the
    # solve.  So a constant sensitivity (e.g. d(R)/d(eta)=0 for additive error) still
    # comes back as a column of that constant -- no special-casing needed.
    .rvarL <- character(0)
    if (!is.null(.rvar)) {
      .rL1 <- c(vapply(dirs[!(dirs %in% .covDirs)], function(.p) paste0("rx_rvar1_", .p, "=", .toRx(.g1(.rvar, .p))), character(1)),
                vapply(.covDirs, function(.cd) .emitScaled1("rx_rvar1_", "", .cd), character(1)))
      .pr2mod <- !(.P2r$i %in% .covDirs) & !(.P2r$j %in% .covDirs)
      .rL2 <- c(vapply(which(.pr2mod), function(.r) .emitPair("rx_rvar2_", .P2r$i[.r], .P2r$j[.r], function(.a, .b) .toRx(.g2(.rvar, .a, .b)), dirs), character(1)),
                vapply(which(!.pr2mod), function(.r) .emitPair("rx_rvar2_", .P2r$i[.r], .P2r$j[.r], NULL, dirs), character(1)))
      .rvarL <- c(paste0("rx_rvarf_=", .toRx(.rvar)), .rL1, .rL2)
    }
    # Residual-variance sigma sensitivities: for each error parameter the variance R
    # depends on (an error-model theta appearing in rx_r_), emit dR/dsigma
    # (rx_rsig_<n>_), its direction sensitivity d2R/(dsigma ddir) (rx_rsig1_<n>_<dir>)
    # and the sigma-pair 2nd derivative d2R/(dsigma dsigma') (rx_rsig2_<n>_<n'>).  In
    # the (f,R) assembly a sigma is a pseudo-direction with df/dsigma=0 and dR/dsigma
    # from here, so every parameter contracts through the same rho(f,R) machinery.
    .sigTh <- integer(0); .sigL <- character(0)
    if (!is.null(.rvar)) {
      .rvarStr <- as.character(.rvar)
      .erTh <- ui$iniDf$ntheta[!is.na(ui$iniDf$err)]
      .frTh <- suppressWarnings(unique(as.integer(sub("THETA_([0-9]+)_", "\\1",
        regmatches(.rvarStr, gregexpr("THETA_[0-9]+_", .rvarStr))[[1]]))))
      .sigTh <- sort(intersect(.erTh, .frTh))
      for (.n in .sigTh) {
        .sg <- paste0("THETA_", .n, "_"); .dRs <- .Dn(.rvar, .sg)
        .sigL <- c(.sigL, paste0("rx_rsig_", .n, "_=", .toRx(.dRs)),
          vapply(dirs[!(dirs %in% .covDirs)], function(.p) paste0("rx_rsig1_", .n, "_", .p, "=", .toRx(.g1(.dRs, .p))), character(1)),
          vapply(.covDirs, function(.cd) .emitScaled1(paste0("rx_rsig1_", .n, "_"), "", .cd), character(1)),
          vapply(.sigTh[.sigTh >= .n], function(.n2)
            paste0("rx_rsig2_", .n, "_", .n2, "=", .toRx(.Dn(.dRs, paste0("THETA_", .n2, "_")))), character(1)))
      }
    }
    # Transform-of-both-sides parameters (rx_yj_/rx_lambda_/rx_low_/rx_hi_): emit copies
    # so the DV can be transformed in R (y = tbs(DV)) the same way the inner model does.
    .tvarL <- character(0); .hasTrans <- FALSE
    for (.tv in c("yj", "lambda", "low", "hi")) {
      .tval <- tryCatch(get(paste0("rx_", .tv, "_"), .s), error = function(e) NULL)
      if (!is.null(.tval)) .tvarL <- c(.tvarL, paste0("rx_t", .tv, "_=", .toRx(.tval)))
    }
    .hasTrans <- length(.tvarL) > 0L && !identical(as.character(ui$predDf$transform), "untransformed")
    .baseOde <- vapply(.st, function(.x) paste0("d/dt(", .x, ")=", .toRx(get(paste0("rx__d_dt_", .x, "__"), .s))), character(1))
    # State initial conditions (state(0) <- expr).  Emit the base IC AND its
    # per-direction 1st/2nd-order sensitivity-compartment ICs -- otherwise the
    # augmented model starts every state (and sensitivity compartment) at 0,
    # which for a model whose prediction is driven by a parameter-dependent IC
    # (e.g. y(0) <- exp(ta), or a DDE history amplitude) makes f -- and every
    # residual-weighted covariance term -- badly wrong.  The IC is a function of
    # parameters/covariates evaluated at t=0 (before integration), so its
    # direction derivatives are direct partials (no state-sensitivity chain).
    # Compartments whose IC .rxSens already emitted (e.g. the DDE delay-sensitivity
    # augmentation writes the sensitivity-compartment histories/ICs itself); skip
    # those to avoid a duplicate `cmt(0)=` assignment.
    .icDone <- trimws(sub("\\(0\\)=.*$", "", grep("\\(0\\)=", unlist(strsplit(c(.s1, .s2), "\n")), value = TRUE)))
    .emitIc <- function(.cmt, .expr) if (identical(.toRx(.expr), "0") || .cmt %in% .icDone) character(0)
      else paste0(.cmt, "(0)=", .toRx(.expr))
    .icL <- character(0)
    for (.x in .st) {
      .ic <- tryCatch(get(paste0("rx_", .x, "_ini_0__"), .s), error = function(e) NULL)
      if (is.null(.ic)) next
      if (!(.x %in% .icDone)) .icL <- c(.icL, paste0(.x, "(0)=", .toRx(.ic)))  # base state IC
      # only model directions have sensitivity compartments (covariate directions are emitted by
      # scaling, so their IC contribution is already in the scaled eta columns -- emitting a
      # phantom rx__sens_..._BY_<covDir>_(0)= for a non-existent compartment would break the model).
      for (.p in .mfDirs) .icL <- c(.icL, .emitIc(paste0("rx__sens_", .x, "_BY_", .p, "__"), .Dn(.ic, .p)))
      .p2m <- which(!(.P2$i %in% .covDirs) & !(.P2$j %in% .covDirs))
      for (.r in .p2m)
        .icL <- c(.icL, .emitIc(paste0("rx__sens_", .x, "_BY_", .P2$i[.r], "_BY_", .P2$j[.r], "__"),
                                .Dn(.Dn(.ic, .P2$i[.r]), .P2$j[.r])))
    }
    # Dosing modifiers (bioavailability f(), lag()/alag(), rate(), dur()) live in the
    # pruned env as rx_<mod>_<state>_ and are NOT part of rx__d_dt_*.  Emit them so the
    # augmented model doses correctly AND so eventSens="jump" (below) fills the analytic
    # dose-parameter ("jump") sensitivities for the rx__sens_* compartments -- otherwise
    # a modeled dosing parameter's sensitivity is silently zero.
    .dosVars <- grep("^rx_(f|lag|alag|rate|dur)_.+_$", ls(envir = .s, all.names = TRUE), value = TRUE)
    .dose <- vapply(.dosVars, function(.v) {
      .m <- regmatches(.v, regexec("^rx_(f|lag|alag|rate|dur)_(.+)_$", .v))[[1]]
      .fun <- if (.m[2L] == "lag") "alag" else .m[2L]     # rxode2 stores lag() as alag()
      paste0(.fun, "(", .m[3L], ")=", .toRx(get(.v, envir = .s)))
    }, character(1))
    # DDE non-constant delay() pre-history: base past(state,tau)<-expr + the
    # per-sensitivity-compartment 1st/2nd-order histories that .rxSens()
    # accumulated as a side effect of the .s1/.s2 calls above (rxode2's
    # .rxDelaySensAugment/.rxDelaySensAugment2).  Spliced after all d/dt +
    # sensitivity-ODE declarations so the referenced state/sens compartments
    # are already defined -- same convention as .rxFinalizeInner (focei.R).
    # NULL (no delay()) for ordinary models.
    .pastLines <- .s$..pastLines
    if (is.null(.pastLines)) .pastLines <- character(0)
    .modTxt <- paste(c(.baseOde, .dose, .s1, .s2, .icL, .pastLines, paste0("rx_predf_=", .toRx(.pred)), .fL1, .fL2, .rvarL, .sigL, .tvarL), collapse = "\n")
    .modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", .modTxt); .modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", .modTxt)
    # Optimize common subexpressions (as the inner model does): the augmented model
    # has heavy shared subexpressions across the sensitivity ODEs and the f1/f2
    # prediction chains, so rxOptExpr materially shrinks the per-solve work.
    # rxode2's rxOptExpr() itself optimizes a large model in cost-balanced chunks by
    # default (it disguises the compartment-scoped constructs -- state ICs and
    # f/alag/lag/rate/dur dosing modifiers -- a chunk would otherwise orphan, and
    # falls back to a single whole-model pass if a chunk cannot be optimized), so the
    # ~200-350 line augmented model no longer needs hand-rolled chunking here.  The
    # fit's rxControl(cores=) is passed through so the chunks are optimized in
    # parallel with the same thread setting the solves use (`.cores` also rides along
    # in the result so the batched solves run parallel).
    .cores <- .optExprCores(ui)
    if (isTRUE(rxode2::rxGetControl(ui, "optExpression", TRUE))) {
      .modTxt <- tryCatch(
        rxode2::rxOptExpr(.modTxt, "FOCEi outer gradient model", parallel = .cores),
        error = function(e) .modTxt)
    }
    # Declare the theta/eta inputs AND the model covariates up front (param()) so the
    # solve parameter order is fixed and positional.  Reuse .uiGetThetaEtaParams -- the
    # SAME theta/eta/covariate ordering the inner sensitivity model uses -- so the
    # augmented outer/analytic model's parameter order matches the inner problem
    # exactly (thetas in ntheta order, etas in neta order, then covariates).
    .param <- .uiGetThetaEtaParams(ui, TRUE)          # params(THETA[1], .., ETA[1], .., covs)
    .param <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", .param))
    .modTxt <- paste(.param, .modTxt, sep = "\n")
    # eventSens="jump" attaches rxode2's analytic event/dosing-parameter sensitivities
    # (forward variational jumps at dose times) for the sensitivity compartments.  `cols`
    # precomputes solve-output column names/index maps; `cores` carries the fit's rxControl thread
    # count so the batched solves run parallel; `key` seeds the per-fit event-table reuse cache.
    list(augMod = rxode2::rxode2(.modTxt, eventSens = "jump"), dirs = dirs, ndir = length(dirs), fDirs = .fDirs,
         st = .st, P2 = .P2, P2r = .P2r, hasRvar = !is.null(.rvar), sigTh = .sigTh, hasTrans = .hasTrans,
         cols = .foceiAnalyticCols(dirs, .fDirs, .P2, .P2r, .sigTh), cores = .cores, key = .key)
  }, error = function(e) NULL)
  if (!is.null(.key) && !is.null(.res)) {
    if (length(ls(.foceiAnalyticAugCache, all.names = TRUE)) >= 64L)    # bound retained compiled models
      rm(list = ls(.foceiAnalyticAugCache, all.names = TRUE), envir = .foceiAnalyticAugCache)
    assign(.key, .res, envir = .foceiAnalyticAugCache)
  }
  .res
}

#' Solve the direction-set 2nd-order model for one subject and recover the
#' 3rd-order tensor `Ath` by Shi (2021) central differences of the analytic
#' 2nd-order sensitivities `A` (C++ `shi21CentralWrap`, perturbing one full-param
#' coordinate per direction), then symmetrize.  Returns `list(f, a, A, Ath)` or
#' `NULL` on failure.  (`rxode2::rxExpandSens3_` would give `Ath` analytically but
#' at O(ndir^3) augmented-model compile cost; the Shi FD keeps it at O(ndir^2).)
#' @noRd
.foceiAnalyticSolveSubjectFD3 <- function(aug, params, ev, times, .fdEps = 7e-7, tol = 1e-10, withR = FALSE) {
  dirs <- aug$dirs; nd <- length(dirs)
  # base solve (f, a, A) via the shared helper; this tier adds the 3rd-order Ath by
  # Shi-differencing A.
  E0 <- .foceiAnalyticSolveFA(aug, params, ev, times, tol = tol); if (is.null(E0)) return(NULL)
  nobs <- length(E0$f)
  # 3rd-order sensitivities Ath = d3f/(deta ddir ddir') and (f,R cov) AthR the same for R,
  # by Shi-differencing the analytic 2nd-order A / AR.  The cov assembly only ever reads
  # these with the FIRST index an ETA (Tn[l,s,t], d2HtDD both have l,m in the eta set), so
  # only the eta directions need to be Shi-differenced -- the theta/sigma "slices" are never
  # used.  This is O(neta) solves instead of O(ndir): the returned tensors are
  # [obs, neta, ndir, ndir] with the eta (differencing) axis first.  A is symmetric in its
  # last two axes, so d(A)/deta is too -- no cross-axis symmetrization is needed.
  # withR=TRUE (the (f,R) cov) also differences AR; the symbolic add/prop cov (R=R(f)) skips it.
  .etaDir <- which(grepl("^ETA_[0-9]+_$", dirs)); neta <- length(.etaDir)
  if (neta == 0L) return(NULL)
  .hasR <- isTRUE(withR) && !is.null(E0$AR); nA <- nobs * nd * nd
  f0 <- if (.hasR) c(as.vector(E0$A), as.vector(E0$AR)) else as.vector(E0$A)
  # shi21CentralWrap differences the closure at a perturbed full-param vector; it
  # strips names, so re-attach names(params) before mapping back through the solve.
  Aflat <- function(.tt) { E <- .foceiAnalyticSolveFA(aug, setNames(.tt, names(params)), ev, times, tol = tol)
    if (is.null(E) || !all(is.finite(E$A))) return(NULL)
    if (.hasR) { if (!all(is.finite(E$AR))) return(NULL); c(as.vector(E$A), as.vector(E$AR)) } else as.vector(E$A) }
  Ath <- array(0, c(nobs, neta, nd, nd)); AthR <- if (.hasR) array(0, c(nobs, neta, nd, nd)) else NULL
  for (li in seq_len(neta)) {
    idx <- match(dirs[.etaDir[li]], names(params))       # this eta's coordinate in params
    if (is.na(idx)) return(NULL)
    sc <- shi21CentralWrap(Aflat, params, f0, idx, .fdEps)  # C++ shi21Central
    if (is.null(sc$gr) || !all(is.finite(sc$gr))) return(NULL)
    Ath[, li, , ] <- array(sc$gr[seq_len(nA)], c(nobs, nd, nd))
    if (.hasR) AthR[, li, , ] <- array(sc$gr[nA + seq_len(nA)], c(nobs, nd, nd))
  }
  .out <- c(list(f = E0$f, a = E0$a, A = E0$A, Ath = Ath),
            if (.hasR) list(R = E0$R, aR = E0$aR, AR = E0$AR, AthR = AthR,
                            Rsig = E0$Rsig, RsigDir = E0$RsigDir, Rsig2 = E0$Rsig2) else NULL,
            if (!is.null(E0$trans)) list(trans = E0$trans) else NULL)
  if (!all(is.finite(.out$f)) || !all(is.finite(.out$a)) || !all(is.finite(.out$A)) || !all(is.finite(.out$Ath))) return(NULL)
  .out
}

#' Route-A (default on; FOCEI_NO_RSIG opts out) sigma reconstruction: the (f,R) cov shares the gradient's `dirs` model
#' (no sigma directions) and rebuilds the sigma tensor slots from the model's own residual-sigma
#' outputs (Rsig=dR/dsigma, RsigDir=d2R/(dsigma ddir), Rsig2=d2R/(dsigma dsigma')) plus their
#' eta-derivatives (RsigDirEta, Rsig2Eta from the FD3).  A sigma has df/dsigma=0, so the prediction
#' slots (a/A/Ath) are zero; the variance slots (aR/AR/AthR) carry the rsig values.  The result is
#' the SAME ndirCov tensor the sigma-as-direction model produced -- but from a model shared with the
#' gradient.
#' @noRd
.foceiAnalyticExpandSigma <- function(E, nsig, neta, RsigDirEta, Rsig2Eta) {
  nd <- ncol(E$a); no <- nrow(E$a); ndc <- nd + nsig; sg <- nd + seq_len(nsig)
  a <- cbind(E$a, matrix(0, no, nsig))
  A <- array(0, c(no, ndc, ndc)); A[, seq_len(nd), seq_len(nd)] <- E$A
  .hasAth <- !is.null(E$Ath)                              # E0 (frozen-R0 solve) has no Ath
  Ath <- if (.hasAth) { .Y <- array(0, c(no, neta, ndc, ndc)); .Y[, , seq_len(nd), seq_len(nd)] <- E$Ath; .Y } else NULL
  aR <- cbind(E$aR, E$Rsig)
  AR <- array(0, c(no, ndc, ndc)); AR[, seq_len(nd), seq_len(nd)] <- E$AR
  .hasAthR <- !is.null(E$AthR) && !is.null(RsigDirEta)     # withR: 3rd-order sigma slots too
  AthR <- if (.hasAthR) { .X <- array(0, c(no, neta, ndc, ndc)); .X[, , seq_len(nd), seq_len(nd)] <- E$AthR; .X } else NULL
  for (k in seq_len(nsig)) {
    AR[, seq_len(nd), sg[k]] <- E$RsigDir[, , k]; AR[, sg[k], seq_len(nd)] <- E$RsigDir[, , k]
    if (.hasAthR) { AthR[, , seq_len(nd), sg[k]] <- RsigDirEta[, , , k]; AthR[, , sg[k], seq_len(nd)] <- RsigDirEta[, , , k] }
    for (l in seq_len(nsig)) { AR[, sg[k], sg[l]] <- E$Rsig2[, k, l]; if (.hasAthR) AthR[, , sg[k], sg[l]] <- Rsig2Eta[, , k, l] }
  }
  E$a <- a; E$A <- A; if (.hasAth) E$Ath <- Ath; E$aR <- aR; E$AR <- AR; if (.hasAthR) E$AthR <- AthR; E
}

#' Batched analogue of [.foceiAnalyticSolveSubjectFD3]: recover the 3rd-order tensor `Ath` (and
#' `AthR` when `withR`) for ALL subjects at once by central-differencing the analytic 2nd-order
#' sensitivities `A` w.r.t. each ETA coordinate via BATCHED population solves ([.foceiAnalyticSolveAll],
#' 1 base + 2*neta perturbed) instead of the per-subject Shi (~O(nsub*neta) solves).  RICHARDSON-
#' extrapolated over steps (h, h/2) to 4th order, with the perturbed solves at a tighter tol, so the
#' batched Ath reproduces the per-subject adaptive-Shi Ath (the FOCEI-vs-FOCE analytic-R agreement,
#' not just the SEs).  Returns the per-subject E-list with `Ath` attached, or NULL (-> per-subject).
#' @noRd
.foceiAnalyticSolveAllFD3 <- function(am, thv, ebes, ids, data, obsTimes, tol = 1e-10,
                                      fdEps = 1e-3, withR = FALSE) {
  dirs <- am$dirs; nd <- length(dirs); neta <- ncol(ebes)
  E0 <- .foceiAnalyticSolveAll(am, thv, ebes, ids, data, obsTimes, tol)
  if (is.null(E0)) return(NULL)
  nsub <- length(E0)
  Ath  <- lapply(E0, function(E) array(0, c(nrow(E$a), neta, nd, nd)))
  AthR <- if (withR) lapply(E0, function(E) array(0, c(nrow(E$a), neta, nd, nd))) else NULL
  # Route A (default ON; FOCEI_NO_RSIG=1 opts out): am is the gradient's `dirs` model (no sigma
  # directions); rebuild the sigma tensor slots from the rsig outputs + their eta-derivatives
  # (differenced here alongside A/AR).  Guarded on Rsig actually being present in the solve.
  .rsig <- !nzchar(Sys.getenv("FOCEI_NO_RSIG")) && !is.null(E0[[1L]]$Rsig) && length(E0[[1L]]$Rsig) > 0L
  nsig <- if (.rsig) ncol(E0[[1L]]$Rsig) else 0L
  RsigDirEta <- if (.rsig && withR) lapply(E0, function(E) array(0, c(nrow(E$a), neta, nd, nsig))) else NULL
  Rsig2Eta   <- if (.rsig && withR) lapply(E0, function(E) array(0, c(nrow(E$a), neta, nsig, nsig))) else NULL
  # differencing A carries the ODE solve's ~tol error floor; solve the PERTURBED models tighter
  # than the base so the 3rd-order Ath stays as clean as the per-subject adaptive Shi.
  .ptol <- min(tol, 1e-12)
  .cdiff <- function(li, h) {
    ep <- ebes; ep[, li] <- ebes[, li] + h; em <- ebes; em[, li] <- ebes[, li] - h
    Ep <- .foceiAnalyticSolveAll(am, thv, ep, ids, data, obsTimes, .ptol)
    Em <- .foceiAnalyticSolveAll(am, thv, em, ids, data, obsTimes, .ptol)
    if (is.null(Ep) || is.null(Em) || length(Ep) != nsub || length(Em) != nsub) return(NULL)
    list(Ep = Ep, Em = Em, h = h)
  }
  for (li in seq_len(neta)) {                             # ETA_li coordinate == ebes column li
    h <- fdEps * max(abs(ebes[, li]), 1)
    s1 <- .cdiff(li, h); s2 <- .cdiff(li, h / 2)
    if (is.null(s1) || is.null(s2)) return(NULL)
    for (i in seq_len(nsub)) {
      d1 <- (s1$Ep[[i]]$A - s1$Em[[i]]$A) / (2 * s1$h); d2 <- (s2$Ep[[i]]$A - s2$Em[[i]]$A) / (2 * s2$h)
      Ath[[i]][, li, , ] <- (4 * d2 - d1) / 3
      if (withR) {
        e1 <- (s1$Ep[[i]]$AR - s1$Em[[i]]$AR) / (2 * s1$h); e2 <- (s2$Ep[[i]]$AR - s2$Em[[i]]$AR) / (2 * s2$h)
        AthR[[i]][, li, , ] <- (4 * e2 - e1) / 3
        if (.rsig) {
          g1 <- (s1$Ep[[i]]$RsigDir - s1$Em[[i]]$RsigDir) / (2 * s1$h); g2 <- (s2$Ep[[i]]$RsigDir - s2$Em[[i]]$RsigDir) / (2 * s2$h)
          RsigDirEta[[i]][, li, , ] <- (4 * g2 - g1) / 3
          k1 <- (s1$Ep[[i]]$Rsig2 - s1$Em[[i]]$Rsig2) / (2 * s1$h); k2 <- (s2$Ep[[i]]$Rsig2 - s2$Em[[i]]$Rsig2) / (2 * s2$h)
          Rsig2Eta[[i]][, li, , ] <- (4 * k2 - k1) / 3
        }
      }
    }
  }
  for (i in seq_len(nsub)) {
    E0[[i]]$Ath <- Ath[[i]]
    if (withR) E0[[i]]$AthR <- AthR[[i]]
    if (.rsig) E0[[i]] <- .foceiAnalyticExpandSigma(E0[[i]], nsig, neta,
                            if (is.null(RsigDirEta)) NULL else RsigDirEta[[i]],
                            if (is.null(Rsig2Eta)) NULL else Rsig2Eta[[i]])
    if (!all(is.finite(E0[[i]]$A)) || !all(is.finite(E0[[i]]$Ath))) return(NULL)
  }
  E0
}

#' Per-subject observed-information R in the general (f,R) form: the prediction f and
#' the variance R are independent solved quantities with sensitivities a/A/Ath and
#' aR/AR/AthR (Ath/AthR by Shi-FD, as the add/prop cov FD-differences A -> Ath).  Every
#' non-Omega parameter (structural theta AND residual sigma) is a DIRECTION; a sigma
#' direction has a=A=Ath=0 (f is sigma-independent) so only its variance sensitivities
#' contribute, and no separate sigma machinery is needed.  The rho(f,R,y) partials up to
#' 3rd order are model-independent closed forms.  `dirP` maps each non-Omega param to its
#' direction; Omega params follow, using `omd`.  Reduces to `.foceiAnalyticSubjectR` when
#' R=R(f).  FOCE (interaction=0) delegates to the frozen-R0 variant.
#' @noRd
.foceiAnalyticSubjectRFR <- function(E, ehat, Om, neta, ndirP, dirP, omd,
                                     ndir = neta, Oi = solve(Om), interaction = 1L,
                                     E0 = NULL, foceType = 0L) {
  if (identical(as.integer(interaction), 0L))
    return(.foceiAnalyticSubjectRfoceFR(E, ehat, Om, neta, ndirP, dirP, omd,
                                        ndir = ndir, Oi = Oi, E0 = E0, foceType = foceType))
  tr <- function(M) sum(diag(M))
  f <- E$f; y <- E$y; R <- E$R; a <- E$a; A <- E$A; Ath <- E$Ath
  aR <- E$aR; AR <- E$AR; AthR <- E$AthR
  res <- y - f
  rf <- -res / R; rR <- 0.5 * (1 / R - res^2 / R^2)
  rff <- 1 / R; rfR <- res / R^2; rRR <- 0.5 * (-1 / R^2 + 2 * res^2 / R^3)
  rffR <- -1 / R^2; rfRR <- -2 * res / R^3; rRRR <- 0.5 * (2 / R^3 - 6 * res^2 / R^4)  # rfff = 0
  iR <- 1 / R; iR2 <- iR^2; iR3 <- iR^3; iR4 <- iR^4
  nom <- omd$nom; np <- ndirP + nom
  ei <- seq_len(neta); di <- seq_len(ndir)
  # (f,R) 2nd total derivative of the density w.r.t. two directions da, db
  Gdd <- function(da, db) sum(rff * a[, da] * a[, db] + rfR * (a[, da] * aR[, db] + aR[, da] * a[, db]) +
                                rRR * aR[, da] * aR[, db] + rf * A[, da, db] + rR * AR[, da, db])
  H <- Oi; for (l in ei) for (m in ei) H[l, m] <- H[l, m] + Gdd(l, m)
  HiM <- solve(H)
  N <- matrix(0, neta, ndir); for (l in ei) for (d in di) N[l, d] <- Gdd(l, d)
  # exact 3rd total derivative Tn[l,s,t] = d2(Phi_l)/ddir_s ddir_t, Phi_l = rf a_l + rR aR_l
  Tn <- array(0, c(neta, ndir, ndir))
  for (l in ei) for (s in di) for (t in di) {
    us  <- rff * a[, s] + rfR * aR[, s]; ut <- rff * a[, t] + rfR * aR[, t]
    ust <- rffR * (a[, s] * aR[, t] + aR[, s] * a[, t]) + rfRR * aR[, s] * aR[, t] + rff * A[, s, t] + rfR * AR[, s, t]
    ws  <- rfR * a[, s] + rRR * aR[, s]; wt <- rfR * a[, t] + rRR * aR[, t]
    wst <- rffR * a[, s] * a[, t] + rfRR * (a[, s] * aR[, t] + aR[, s] * a[, t]) + rRRR * aR[, s] * aR[, t] +
      rfR * A[, s, t] + rRR * AR[, s, t]
    Tn[l, s, t] <- sum(ust * a[, l] + us * A[, l, t] + ut * A[, l, s] + rf * Ath[, l, s, t] +
                         wst * aR[, l] + ws * AR[, l, t] + wt * AR[, l, s] + rR * AthR[, l, s, t])
  }
  # Laplace determinant Ht = Oi + sum(a a / R + 0.5 aR aR / R^2) and its 1st/2nd
  # direction derivatives (the interaction-free E[rho_fR]=0 split into f- and R-quadratics)
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(a[, l] * a[, m] * iR + 0.5 * aR[, l] * aR[, m] * iR2)
  Hti <- solve(Ht)
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum((A[, l, s] * a[, m] + a[, l] * A[, m, s]) * iR - a[, l] * a[, m] * aR[, s] * iR2 +
                     0.5 * (AR[, l, s] * aR[, m] + aR[, l] * AR[, m, s]) * iR2 - aR[, l] * aR[, m] * aR[, s] * iR3); D })
  d2HtDD <- lapply(di, function(s) lapply(di, function(t) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(
      (Ath[, l, s, t] * a[, m] + A[, l, s] * A[, m, t] + A[, l, t] * A[, m, s] + a[, l] * Ath[, m, s, t]) * iR -
        (A[, l, s] * a[, m] + a[, l] * A[, m, s]) * aR[, t] * iR2 -
        (A[, l, t] * a[, m] + a[, l] * A[, m, t]) * aR[, s] * iR2 - a[, l] * a[, m] * AR[, s, t] * iR2 +
        2 * a[, l] * a[, m] * aR[, s] * aR[, t] * iR3 +
        0.5 * (AthR[, l, s, t] * aR[, m] + AR[, l, s] * AR[, m, t] + AR[, l, t] * AR[, m, s] + aR[, l] * AthR[, m, s, t]) * iR2 -
        (AR[, l, s] * aR[, m] + aR[, l] * AR[, m, s]) * aR[, t] * iR3 -
        (AR[, l, t] * aR[, m] + aR[, l] * AR[, m, t]) * aR[, s] * iR3 - aR[, l] * aR[, m] * AR[, s, t] * iR3 +
        3 * aR[, l] * aR[, m] * aR[, s] * aR[, t] * iR4); D }))
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtD[[l]]), numeric(1))
  Cee <- matrix(0, neta, neta); for (s in ei) for (t in ei)
    Cee[s, t] <- 0.5 * (tr(Hti %*% d2HtDD[[s]][[t]]) - tr(Hti %*% dHtD[[s]] %*% Hti %*% dHtD[[t]]))
  typ <- function(p) if (p <= ndirP) "dir" else "om"
  dOf <- function(p) dirP[p]; omc <- function(p) p - ndirP
  Mcol <- function(p) if (typ(p) == "dir") N[, dOf(p)] else as.numeric(omd$dOi[[omc(p)]] %*% ehat)
  dHt_p <- function(p) if (typ(p) == "dir") dHtD[[dOf(p)]] else omd$dOi[[omc(p)]]
  d2HtEtaP <- function(p, l) if (typ(p) == "dir") d2HtDD[[dOf(p)]][[l]] else matrix(0, neta, neta)
  d2Ht_pp <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "dir" && tb == "dir") return(d2HtDD[[dOf(aa)]][[dOf(bb)]])
    if (ta == "om" && tb == "om") return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    matrix(0, neta, neta) }
  Smat <- function(p) { if (typ(p) == "om") return(omd$dOi[[omc(p)]])
    M <- matrix(0, neta, ndir); for (l in ei) for (s in di) M[l, s] <- Tn[l, dOf(p), s]; M }
  Svec <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "dir" && tb == "dir") return(Tn[, dOf(aa), dOf(bb)])
    if (ta == "om" && tb == "om") return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    rep(0, neta) }
  d2Phi <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "dir" && tb == "dir") return(Gdd(dOf(aa), dOf(bb)))
    if (ta == "om" && tb == "om") return(0.5 * as.numeric(t(ehat) %*% omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat) + 0.5 * omd$d2LD[omc(aa), omc(bb)])
    0 }
  .Mcols <- lapply(1:np, Mcol)
  etaP <- matrix(vapply(1:np, function(p) as.numeric(-HiM %*% .Mcols[[p]]), numeric(neta)), nrow = neta)
  eta2 <- function(aa, bb) { b <- Svec(aa, bb) + Smat(aa)[, ei, drop = FALSE] %*% etaP[, bb] + Smat(bb)[, ei, drop = FALSE] %*% etaP[, aa]
    for (l in ei) b[l] <- b[l] + as.numeric(t(etaP[, aa]) %*% Tn[l, ei, ei] %*% etaP[, bb]); as.numeric(-HiM %*% b) }
  Cpe <- function(p, l) 0.5 * (tr(Hti %*% d2HtEtaP(p, l)) - tr(Hti %*% dHt_p(p) %*% Hti %*% dHtD[[l]]))
  Cpp <- function(aa, bb) 0.5 * (tr(Hti %*% d2Ht_pp(aa, bb)) - tr(Hti %*% dHt_p(aa) %*% Hti %*% dHt_p(bb)))
  .CpeRow <- lapply(1:np, function(p) vapply(ei, function(l) Cpe(p, l), numeric(1)))
  R <- matrix(0, np, np)
  for (aa in 1:np) for (bb in aa:np) {
    dat <- d2Phi(aa, bb) - as.numeric(t(.Mcols[[aa]]) %*% HiM %*% .Mcols[[bb]])
    ld <- Cpp(aa, bb) + sum(.CpeRow[[aa]] * etaP[, bb]) + sum(.CpeRow[[bb]] * etaP[, aa]) +
          as.numeric(t(etaP[, aa]) %*% Cee %*% etaP[, bb]) + sum(Cen * eta2(aa, bb))
    R[aa, bb] <- R[bb, aa] <- dat + ld
  }
  R
}

#' (f,R) FOCE per-subject observed-information R (interaction=0).  The inner problem is
#' interaction-free (S_FOCE = sum(q0 a) + Omega^-1 eta, q0=-(y-f)/R0, q1=1/R0), so the
#' EBE does not stationarize the full Laplace objective and the general non-envelope form
#' R_ab = F_ab + F_aeta eta_b + F_beta eta_a + eta_a' F_etaeta eta_b + F_eta eta_ab is used.
#' R0 is frozen at the eta=0 population variance (nonmem, from E0) or the live conditional
#' variance (foce+, from E); its theta-chain enters the parameter columns via aRc/ARc
#' (dR0/ddir, d2R0/ddir2 from E0), with the eta-block frozen (aRe=0 for nonmem).  Ath is
#' the eta-hat prediction 3rd-order tensor [obs, neta, ndir, ndir].  Sigmas are directions
#' (a=A=Ath=0).  Reduces to `.foceiAnalyticSubjectRfoce` on add/prop.
#' @noRd
.foceiAnalyticSubjectRfoceFR <- function(E, ehat, Om, neta, ndirP, dirP, omd,
                                         ndir = neta, Oi = solve(Om), E0 = NULL, foceType = 0L) {
  tr <- function(M) sum(diag(M))
  f <- E$f; y <- E$y; a <- E$a; A <- E$A; Ath <- E$Ath
  ei <- seq_len(neta); di <- seq_len(ndir); nobs <- length(f)
  .fp <- identical(as.integer(foceType), 1L) || is.null(E0)
  # frozen variance R0 and its sensitivities.  aRe drives the eta-block (0 for nonmem --
  # R0 is eta-independent -- live E$aR for foce+); aRc drives the parameter columns
  # (E0's full dR0/ddir for nonmem, live E$aR for foce+), so a mu-referenced theta (whose
  # direction is an eta) keeps its frozen-R0 theta sensitivity.
  if (.fp) { R0 <- E$R; aRe <- E$aR; aRc <- E$aR; ARc <- E$AR }
  else { R0 <- E0$R; aRe <- matrix(0, nobs, ndir); aRc <- E0$aR; ARc <- E0$AR }
  res <- y - f
  rf <- -res / R0; rR <- 0.5 * (1 / R0 - res^2 / R0^2)
  rff <- 1 / R0; rfR <- res / R0^2; rRR <- 0.5 * (-1 / R0^2 + 2 * res^2 / R0^3)
  q0 <- rf; q1 <- 1 / R0; iR <- 1 / R0; iR2 <- iR^2; iR3 <- iR^3
  nom <- omd$nom; np <- ndirP + nom
  ae <- a[, ei, drop = FALSE]
  isD <- function(p) p <= ndirP; dOf <- function(p) dirP[p]; omc <- function(p) p - ndirP
  # ---- Phi (data) tensors: H=Phi_etaeta, gPhi=Phi_eta (aRe eta-block) ----
  gPhi <- as.numeric(Oi %*% ehat); for (l in ei) gPhi[l] <- gPhi[l] + sum(rf * a[, l] + rR * aRe[, l])
  H <- Oi; for (l in ei) for (m in ei)
    H[l, m] <- H[l, m] + sum(rff * a[, l] * a[, m] + rfR * (a[, l] * aRe[, m] + aRe[, l] * a[, m]) +
                              rRR * aRe[, l] * aRe[, m] + rf * A[, l, m] + rR * E_ARelm(E, l, m, .fp))
  # ---- FOCE inner (EBE) tensors: interaction-free q-based Hf/Nf/Tnf ----
  Hf <- Oi; Nf <- matrix(0, neta, ndir)
  for (l in ei) { for (m in ei) Hf[l, m] <- Hf[l, m] + sum(q1 * a[, l] * a[, m] + q0 * A[, l, m])
    for (d in di) Nf[l, d] <- sum(q1 * a[, l] * a[, d] + q0 * A[, l, d]) }
  HfInv <- solve(Hf)
  Tnf <- array(0, c(neta, ndir, ndir)); for (l in ei) for (s in di) for (t in di)
    Tnf[l, s, t] <- sum(q1 * (A[, l, s] * a[, t] + A[, l, t] * a[, s] + A[, s, t] * a[, l]) + q0 * Ath[, l, s, t])
  # ---- determinant Ht = Oi + sum(a a / R0) (interaction-free) + its derivatives ----
  # dHtDir/d2HtDir are parameterized by the R0-sensitivity of each direction: the eta-block
  # uses aRe (frozen, 0 for nonmem), the parameter columns use aRc (E0's dR0/ddir), so a
  # mu-referenced theta keeps its frozen-R0 chain.  ARv is d2R0 for the pair (0 unless both
  # indices carry a live R0 dependence).
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(a[, l] * a[, m] * iR); Hti <- solve(Ht)
  ARblk <- function(s, t) if (.fp) E$AR[, s, t] else rep(0, nobs)   # d2R0/(dir_s dir_t), eta-block
  dHtDir <- function(s, aRvS) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum((A[, l, s] * a[, m] + a[, l] * A[, m, s]) * iR - a[, l] * a[, m] * aRvS[, s] * iR2); D }
  d2HtDir <- function(s, t, aRvS, aRvT, ARv) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum((Ath[, l, s, t] * a[, m] + A[, l, s] * A[, m, t] + A[, l, t] * A[, m, s] + a[, l] * Ath[, m, s, t]) * iR -
      (A[, l, s] * a[, m] + a[, l] * A[, m, s]) * aRvT[, t] * iR2 -
      (A[, l, t] * a[, m] + a[, l] * A[, m, t]) * aRvS[, s] * iR2 - a[, l] * a[, m] * ARv * iR2 +
      2 * a[, l] * a[, m] * aRvS[, s] * aRvT[, t] * iR3); D }
  dHtE <- lapply(ei, function(l) dHtDir(l, aRe))                    # eta-block dHt/deta_l
  d2HtEE <- lapply(ei, function(s) lapply(ei, function(t) d2HtDir(s, t, aRe, aRe, ARblk(s, t))))
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtE[[l]]), numeric(1))
  Cee <- matrix(0, neta, neta); for (s in ei) for (t in ei)
    Cee[s, t] <- 0.5 * (tr(Hti %*% d2HtEE[[s]][[t]]) - tr(Hti %*% dHtE[[s]] %*% Hti %*% dHtE[[t]]))
  dHtP <- function(p) if (isD(p)) dHtDir(dOf(p), aRc) else omd$dOi[[omc(p)]]  # parameter dHt/dp
  # ---- R0 theta-chains for the parameter accessors (aRc/ARc); sigma is a direction ----
  chQ <- function(d) (res / R0^2) * aRc[, d]                       # d(q0)/dtheta via R0
  # Phi_(eta,p) and S_p share the (res/R0^2) aRc chain (q0 = Phi_f for FOCE)
  McolData <- function(p) { if (!isD(p)) return(as.numeric(omd$dOi[[omc(p)]] %*% ehat))
    d <- dOf(p); Nf[, d] + as.numeric(crossprod(ae, chQ(d))) }
  McolEBE <- McolData
  d2Phi <- function(aa, bb) { if (!isD(aa) && !isD(bb))
      return(0.5 * as.numeric(t(ehat) %*% omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat) + 0.5 * omd$d2LD[omc(aa), omc(bb)])
    if (!isD(aa) || !isD(bb)) return(0)
    da <- dOf(aa); db <- dOf(bb)
    sum(rff * a[, da] * a[, db] + rf * A[, da, db] + rfR * (a[, da] * aRc[, db] + aRc[, da] * a[, db]) +
          rRR * aRc[, da] * aRc[, db] + rR * ARc[, da, db]) }
  # S_(p,eta) row (SmatEBE): Tnf plus the theta chain -aRc/R0^2 a a + (res/R0^2) aRc A
  SmatEBE <- function(p) { if (!isD(p)) return(omd$dOi[[omc(p)]])
    d <- dOf(p); M <- matrix(0, neta, ndir); for (l in ei) for (s in di)
      M[l, s] <- Tnf[l, d, s] + sum(-aRc[, d] * iR2 * a[, s] * a[, l] + (res * iR2) * aRc[, d] * A[, l, s]); M }
  # S_(p,p') vector (SvecEBE): Tnf + the combined 2nd-order R0 chain (R0'A0 cancels)
  SvecEBE <- function(aa, bb) { ta <- isD(aa); tb <- isD(bb)
    if (!ta && !tb) return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    if (!ta || !tb) return(rep(0, neta))
    da <- dOf(aa); db <- dOf(bb); v <- Tnf[, da, db]
    w <- -(a[, da] * aRc[, db] + aRc[, da] * a[, db]) * iR2 + (res * iR2) * ARc[, da, db] - 2 * res * aRc[, da] * aRc[, db] * iR3
    for (l in ei) v[l] <- v[l] + sum(w * a[, l] + (res * iR2) * (aRc[, da] * A[, l, db] + aRc[, db] * A[, l, da])); v }
  # determinant d2Ht/(deta_l dp) uses aRc for the theta-direction, aRe for the eta; the mixed
  # d2R0/(deta dtheta) is 0 for nonmem (R0 frozen w.r.t. eta) and E$AR for foce+.
  d2HtEtaP <- function(p, l) { if (!isD(p)) return(matrix(0, neta, neta))
    d <- dOf(p); d2HtDir(d, l, aRc, aRe, if (.fp) E$AR[, d, l] else rep(0, nobs)) }
  d2Ht_pp <- function(aa, bb) { ta <- isD(aa); tb <- isD(bb)
    if (ta && tb) return(d2HtDir(dOf(aa), dOf(bb), aRc, aRc, ARc[, dOf(aa), dOf(bb)]))
    if (!ta && !tb) return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    matrix(0, neta, neta) }
  Cpe <- function(p, l) 0.5 * (tr(Hti %*% d2HtEtaP(p, l)) - tr(Hti %*% dHtP(p) %*% Hti %*% dHtE[[l]]))
  Cpp <- function(aa, bb) 0.5 * (tr(Hti %*% d2Ht_pp(aa, bb)) -
    tr(Hti %*% dHtP(aa) %*% Hti %*% dHtP(bb)))
  .McD <- lapply(1:np, McolData)
  etaP <- matrix(vapply(1:np, function(p) as.numeric(-HfInv %*% McolEBE(p)), numeric(neta)), nrow = neta)
  eta2 <- function(aa, bb) { b <- SvecEBE(aa, bb) + SmatEBE(aa)[, ei, drop = FALSE] %*% etaP[, bb] +
      SmatEBE(bb)[, ei, drop = FALSE] %*% etaP[, aa]
    for (l in ei) b[l] <- b[l] + as.numeric(t(etaP[, aa]) %*% Tnf[l, ei, ei] %*% etaP[, bb]); as.numeric(-HfInv %*% b) }
  .CpeRow <- lapply(1:np, function(p) vapply(ei, function(l) Cpe(p, l), numeric(1)))
  R <- matrix(0, np, np)
  for (aa in 1:np) for (bb in aa:np) {
    e_ab <- eta2(aa, bb)
    dat <- d2Phi(aa, bb) + sum(.McD[[aa]] * etaP[, bb]) + sum(.McD[[bb]] * etaP[, aa]) +
           as.numeric(t(etaP[, aa]) %*% H %*% etaP[, bb]) + sum(gPhi * e_ab)
    ld <- Cpp(aa, bb) + sum(.CpeRow[[aa]] * etaP[, bb]) + sum(.CpeRow[[bb]] * etaP[, aa]) +
          as.numeric(t(etaP[, aa]) %*% Cee %*% etaP[, bb]) + sum(Cen * e_ab)
    R[aa, bb] <- R[bb, aa] <- dat + ld
  }
  R
}

#' d2R0/(deta_l deta_m) for the FOCE eta-block H: 0 (nonmem, R0 frozen) or live E$AR (foce+).
#' @noRd
E_ARelm <- function(E, l, m, fp) if (fp) E$AR[, l, m] else 0

#' C++/Armadillo port of `.foceiAnalyticSubjectRfoceFR` (FOCE (f,R) observed information).
#' Resolves the frozen-R0 sensitivities (eta-block aRe/ARe: 0 for nonmem, live E for foce+;
#' parameter columns aRc/ARc: E0 for nonmem, live E for foce+), reshapes Ath, and calls the
#' kernel.  Matches `.foceiAnalyticSubjectRfoceFR` exactly.
#' @noRd
.foceiAnalyticSubjectRfoceFRCpp <- function(E, ehat, Om, neta, ndirP, dirP, omd, ndir,
                                            Oi = solve(Om), E0 = NULL, foceType = 0L,
                                            dvSens = matrix(0, length(E$f), 0L),
                                            dvSens2 = matrix(0, length(E$f), 0L),
                                            censv = integer(0), limv = numeric(0)) {
  nobs <- length(E$f); nom <- omd$nom
  .fp <- identical(as.integer(foceType), 1L) || is.null(E0)
  if (.fp) { R0 <- E$R; aRe <- E$aR; aRc <- E$aR; ARc <- E$AR; ARe <- E$AR }
  else { R0 <- E0$R; aRe <- matrix(0, nobs, ndir); aRc <- E0$aR; ARc <- E0$AR; ARe <- array(0, c(nobs, ndir, ndir)) }
  AthC <- array(E$Ath, c(nobs, neta, ndir * ndir))
  dOiC <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiC[, , k] <- omd$dOi[[k]]
  d2OiC <- array(0, c(neta, neta, max(nom * nom, 1L)))
  if (nom > 0L) for (aa in seq_len(nom)) for (bb in seq_len(nom)) d2OiC[, , (aa - 1L) * nom + bb] <- omd$d2Oi[[aa]][[bb]]
  d2LD <- if (nom > 0L) omd$d2LD else matrix(0, 1, 1)
  foceiSubjectRfoceFR_(E$a, E$A, AthC, aRe, aRc, ARe, ARc, dvSens, dvSens2, as.integer(censv), as.numeric(limv),
                       E$f, E$y, R0, as.numeric(ehat), Oi,
                       dOiC, d2OiC, d2LD, neta, ndir, ndirP, nom, as.integer(dirP))
}

#' C++/Armadillo port of `.foceiAnalyticSubjectRFR` (FOCEI (f,R) observed information).
#' Reshapes the 3rd-order Ath/AthR tensors and the Omega derivative lists for the kernel.
#' @noRd
.foceiAnalyticSubjectRFRCpp <- function(E, ehat, Om, neta, ndirP, dirP, omd, ndir, Oi = solve(Om),
                                        dvSens = matrix(0, length(E$f), 0L),
                                        dvSens2 = matrix(0, length(E$f), 0L),
                                        censv = integer(0), limv = numeric(0)) {
  nobs <- length(E$f); nom <- omd$nom
  # Ath/AthR are [obs, neta, ndir, ndir] (eta axis first); reshape to (obs, neta, ndir^2)
  # so the kernel reads Ath(o, l, s + t*ndir) with l over the etas.
  AthC <- array(E$Ath, c(nobs, neta, ndir * ndir)); AthRC <- array(E$AthR, c(nobs, neta, ndir * ndir))
  dOiC <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiC[, , k] <- omd$dOi[[k]]
  d2OiC <- array(0, c(neta, neta, max(nom * nom, 1L)))
  if (nom > 0L) for (aa in seq_len(nom)) for (bb in seq_len(nom)) d2OiC[, , (aa - 1L) * nom + bb] <- omd$d2Oi[[aa]][[bb]]
  d2LD <- if (nom > 0L) omd$d2LD else matrix(0, 1, 1)
  foceiSubjectRFR_(E$a, E$A, AthC, E$aR, E$AR, AthRC, dvSens, dvSens2, as.integer(censv), as.numeric(limv),
                   E$f, E$y, E$R, as.numeric(ehat), Oi,
                   dOiC, d2OiC, d2LD, neta, ndir, ndirP, nom, as.integer(dirP))
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
                                   interaction = 1L, E0 = NULL, foceType = 0L) {
  if (identical(as.integer(interaction), 0L))
    return(.foceiAnalyticSubjectRfoce(E, ehat, Om, ef, neta, nth, nsg, sgVar, omd,
                                      ndir = ndir, dirTh = dirTh, Oi = Oi, E0 = E0,
                                      foceType = foceType))
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
                                       ndir = neta, dirTh = seq_len(nth), Oi = solve(Om), E0 = NULL,
                                       foceType = 0L) {
  tr <- function(M) sum(diag(M))
  a <- E$a; A <- E$A; Ath <- E$Ath; f <- E$f; y <- E$y
  # CORRECTED FOCE (foceType=0): the variance R0 is frozen at the eta=0 POPULATION
  # prediction f0 (E0$f), with population sensitivities a0/A0 (E0$a/E0$A) supplying R0's
  # theta-chain.  The residual/numerator/inner sensitivities (y-f, f, a, A, Ath) stay at
  # eta-hat.  R0 is eta-independent (dR0/deta=0), so q1=1/R0, q2=0, pF1=0, r3=0; its
  # theta/sigma derivatives (through f0) enter as a0-chain corrections on the `th`
  # accessors only.  "foce+" (foceType=1) keeps the live conditional R: the live-R
  # pieces carry the dR/df chain through f's own sensitivities, so no f0 solve/chain.
  .fc <- if (identical(as.integer(foceType), 1L)) ef$focePlus else ef$foce
  f0 <- if (!is.null(E0)) E0$f else f
  evf <- function(e) ef$ev(e, f, y, f0)
  rd <- list(r1 = evf(.fc$sc$r1), r2 = evf(.fc$sc$r2), r3 = evf(.fc$sc$r3))   # full rho (Phi), R0
  qd <- list(q0 = evf(.fc$sc$q0), q1 = evf(.fc$sc$q1), q2 = evf(.fc$sc$q2))   # FOCE inner gradient coef, R0
  pf <- list(p = evf(.fc$sc$pF), p1 = evf(.fc$sc$pF1), p2 = evf(.fc$sc$pF2))  # FOCE determinant p = 1/R0
  np <- nth + nsg + omd$nom
  ei <- seq_len(neta); di <- seq_len(ndir); ae <- a[, ei, drop = FALSE]
  .dirOf <- function(p) dirTh[p]
  # a0-chain (only when the variance depends on the population prediction, i.e. prop part)
  .cf0 <- isTRUE(.fc$dependsF0) && !is.null(E0)
  a0 <- if (.cf0) E0$a else NULL; A0 <- if (.cf0) E0$A else NULL
  fq <- if (.cf0) list(qf0 = evf(.fc$f0$qf0), qff0 = evf(.fc$f0$qff0), qf0f0 = evf(.fc$f0$qf0f0),
                       pFf0 = evf(.fc$f0$pFf0), pFf0f0 = evf(.fc$f0$pFf0f0),
                       rhof0 = evf(.fc$f0$rhof0), rhof0f0 = evf(.fc$f0$rhof0f0)) else NULL

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
  perR <- function(p) lapply(.fc$per[[sgi(p)]][c("rf", "rff")], evf)                 # rho sigma partials (Phi), R0
  perQ <- function(p) list(qs = evf(.fc$per[[sgi(p)]]$qs), qsf = evf(.fc$per[[sgi(p)]]$qsf))  # FOCE inner sigma, R0
  perP <- function(p) list(ps = evf(.fc$per[[sgi(p)]]$psF), pf = evf(.fc$per[[sgi(p)]]$pfF))   # FOCE det sigma, R0
  # sigma x population-prediction cross fields (a0-chain of the sigma x theta blocks)
  PVf0 <- function(p) lapply(.fc$perf0[[sgi(p)]], evf)                               # qf0s, pFf0s, rhof0s
  pairKey <- function(aa, bb) { s1 <- sgi(aa); s2 <- sgi(bb); if (paste0(s1, s2) %in% names(.fc$pair)) paste0(s1, s2) else paste0(s2, s1) }
  pairR <- function(aa, bb) lapply(.fc$pair[[pairKey(aa, bb)]][c("rss", "rfss")], evf)
  pairQ <- function(aa, bb) list(qss = evf(.fc$pair[[pairKey(aa, bb)]]$qss))
  pairP <- function(aa, bb) list(pss = evf(.fc$pair[[pairKey(aa, bb)]]$pssF))

  # ---- data-side accessors (Phi partials) ----
  # a0-correction to a theta-derivative of q0=Phi_f: d(q0)/d(theta) picks up qf0*a0 (the
  # population prediction moves with theta); rho_ff0 == qf0.
  McolData <- function(p) { t <- typ(p)                                            # Phi_(eta,p)
    if (t == "th") { v <- Ndat[, .dirOf(p)]
      if (.cf0) v <- v + as.numeric(crossprod(ae, fq$qf0 * a0[, .dirOf(p)])); return(v) }
    if (t == "sg") return(as.numeric(crossprod(ae, perR(p)$rf)))
    as.numeric(omd$dOi[[omc(p)]] %*% ehat) }
  d2Phi <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)                          # Phi_(p,p) at fixed eta
    if (ta == "th" && tb == "th") { da <- .dirOf(aa); db <- .dirOf(bb)
      v <- sum(rd$r2 * a[, da] * a[, db] + rd$r1 * A[, da, db])
      if (.cf0) v <- v + sum(fq$qf0 * (a[, da] * a0[, db] + a0[, da] * a[, db]) +
                             fq$rhof0f0 * a0[, da] * a0[, db] + fq$rhof0 * A0[, da, db])
      return(v) }
    if (ta == "om" && tb == "om") return(0.5 * as.numeric(t(ehat) %*% omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat) + 0.5 * omd$d2LD[omc(aa), omc(bb)])
    if (ta == "om" || tb == "om") return(0)
    if (ta == "sg" && tb == "sg") return(sum(pairR(aa, bb)$rss))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa
    v <- as.numeric(crossprod(a[, .dirOf(thp)], perR(sg)$rf))
    if (.cf0) v <- v + as.numeric(crossprod(a0[, .dirOf(thp)], PVf0(sg)$rhof0s)); v }

  # ---- EBE-side accessors (S_FOCE partials); S=Phi_eta, so S_p carries the same qf0*a0
  # theta-chain as McolData; the eta-role indices (contracted with etaP) get no a0 term ----
  McolEBE <- function(p) { t <- typ(p)                                             # S_p
    if (t == "th") { v <- Nf[, .dirOf(p)]
      if (.cf0) v <- v + as.numeric(crossprod(ae, fq$qf0 * a0[, .dirOf(p)])); return(v) }
    if (t == "sg") return(as.numeric(crossprod(ae, perQ(p)$qs)))
    as.numeric(omd$dOi[[omc(p)]] %*% ehat) }
  SmatEBE <- function(p) { t <- typ(p)                                             # S_(p,eta), eta-role s
    if (t == "th") { d <- .dirOf(p); M <- matrix(0, neta, ndir)
      for (l in ei) for (s in di) M[l, s] <- Tnf[l, d, s]
      if (.cf0) for (l in ei) for (s in ei)
        M[l, s] <- M[l, s] + sum(fq$qff0 * a[, s] * a0[, d] * a[, l] + fq$qf0 * a0[, d] * A[, l, s])
      return(M) }
    if (t == "om") return(omd$dOi[[omc(p)]])
    P <- perQ(p); M <- matrix(0, neta, ndir); for (l in ei) for (s in di) M[l, s] <- sum(P$qsf * a[, s] * a[, l] + P$qs * A[, l, s]); M }
  SvecEBE <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)                       # S_(p,p)
    if (ta == "th" && tb == "th") { da <- .dirOf(aa); db <- .dirOf(bb); v <- Tnf[, da, db]
      if (.cf0) { w <- fq$qff0 * (a[, da] * a0[, db] + a0[, da] * a[, db]) + fq$qf0f0 * a0[, da] * a0[, db] + fq$qf0 * A0[, da, db]
        for (l in ei) v[l] <- v[l] + sum(w * a[, l] + fq$qf0 * (a0[, da] * A[, l, db] + a0[, db] * A[, l, da])) }
      return(v) }
    if (ta == "om" && tb == "om") return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    if (ta == "om" || tb == "om") return(rep(0, neta))
    if (ta == "sg" && tb == "sg") return(as.numeric(crossprod(ae, pairQ(aa, bb)$qss)))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa
    v <- SmatEBE(sg)[, .dirOf(thp)]
    if (.cf0) v <- v + as.numeric(crossprod(ae, PVf0(sg)$qf0s * a0[, .dirOf(thp)])); v }

  # ---- determinant-side accessors (H~ partials, p = 1/R0).  dHtD/d2HtDD are eta-role
  # (pF1=pF2=0), so a theta-derivative of pF=1/R0 adds a pFf0*a0 chain on the `th` blocks. ----
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") { D <- dHtD[[.dirOf(p)]]
      if (.cf0) { d <- .dirOf(p); for (m in ei) for (n in ei) D[m, n] <- D[m, n] + sum(fq$pFf0 * a0[, d] * a[, m] * a[, n]) }
      return(D) }
    if (t == "sg") return(ouAA(perP(p)$ps)); omd$dOi[[omc(p)]] }
  d2HtEtaP <- function(p, l) { t <- typ(p)
    if (t == "th") { D <- d2HtDD[[.dirOf(p)]][[l]]
      if (.cf0) { d <- .dirOf(p); for (m in ei) for (n in ei)
        D[m, n] <- D[m, n] + sum(fq$pFf0 * a0[, d] * (A[, m, l] * a[, n] + a[, m] * A[, n, l])) }
      return(D) }
    if (t == "om") return(matrix(0, neta, neta))
    P <- perP(p); D <- matrix(0, neta, neta); for (s in ei) for (m in ei)
      D[s, m] <- sum(P$pf * a[, l] * a[, s] * a[, m] + P$ps * (A[, s, l] * a[, m] + a[, s] * A[, m, l])); D }
  d2Ht_pp <- function(aa, bb) { ta <- typ(aa); tb <- typ(bb)
    if (ta == "th" && tb == "th") { da <- .dirOf(aa); db <- .dirOf(bb); D <- d2HtDD[[da]][[db]]
      if (.cf0) for (m in ei) for (n in ei)
        D[m, n] <- D[m, n] + sum(fq$pFf0f0 * a0[, da] * a0[, db] * a[, m] * a[, n] + fq$pFf0 * A0[, da, db] * a[, m] * a[, n] +
          fq$pFf0 * (a0[, da] * A[, m, db] + a0[, db] * A[, m, da]) * a[, n] +
          fq$pFf0 * (a0[, da] * A[, n, db] + a0[, db] * A[, n, da]) * a[, m])
      return(D) }
    if (ta == "om" && tb == "om") return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    if (ta == "om" || tb == "om") return(matrix(0, neta, neta))
    if (ta == "sg" && tb == "sg") return(ouAA(pairP(aa, bb)$pss))
    thp <- if (ta == "th") aa else bb; sg <- if (ta == "th") bb else aa; d <- .dirOf(thp)
    D <- d2HtEtaP(sg, d)
    if (.cf0) { P <- PVf0(sg); for (m in ei) for (n in ei) D[m, n] <- D[m, n] + sum(P$pFf0s * a0[, d] * a[, m] * a[, n]) }
    D }
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
  .ev <- .foceiAnalyticEvents(aug, ev)                    # reuse the pre-translated events (FOCE Newton reuse)
  .nc <- if (is.null(aug$cores)) 0L else aug$cores
  # DDE augmented solve: force pure dop853 (dense, no Jacobian).  Its 8th-order dense history
  # reproduces the delayed sensitivity solve exactly and needs no Jacobian, so it sidesteps the
  # composite/ros4 on-the-fly Jacobian generation for this THETA/ETA-named augmented model.
  .ddeArgs <- if (isTRUE(rxode2::rxModelVars(aug$augMod)$flags[["hasDelay"]] == 1L))
    list(method = "dop853", stiff2 = 0L, dense = TRUE) else list()
  .d <- tryCatch(withCallingHandlers(
      as.data.frame(do.call(rxode2::rxSolve, c(list(aug$augMod, params = params, events = .ev, cores = .nc,
          returnType = "data.frame", atol = tol, rtol = tol), .ddeArgs))),
      warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) NULL)
  if (is.null(.d)) return(NULL)
  .d <- .d[.d$time %in% times, , drop = FALSE]
  .fD <- if (is.null(aug$fDirs)) dirs else aug$fDirs   # sigma-skip: f-sensitivities only over f-directions
  if (nrow(.d) == 0L || nrow(.d) != length(times) ||
        !all(c("rx_predf_", paste0("rx_f1_", .fD)) %in% names(.d))) return(NULL)
  no <- nrow(.d)
  a <- matrix(0, no, nd); a[, match(.fD, dirs)] <- vapply(.fD, function(q) .d[[paste0("rx_f1_", q)]], numeric(no))
  A <- array(0, c(no, nd, nd))
  for (r in seq_len(nrow(aug$P2))) {                   # only i<=j emitted -> mirror to i>j (f-dir pairs)
    .ii <- match(aug$P2$i[r], dirs); .jj <- match(aug$P2$j[r], dirs)
    .v2 <- .d[[paste0("rx_f2_", aug$P2$i[r], "_", aug$P2$j[r])]]
    A[, .ii, .jj] <- .v2; A[, .jj, .ii] <- .v2
  }
  .out <- list(f = .d$rx_predf_, a = a, A = A)
  if (isTRUE(aug$hasRvar)) .out <- c(.out, .foceiReadRvar(.d, aug, no))
  .out
}

#' Read the residual variance R and its 1st/2nd direction sensitivities (aR, AR)
#' from an augmented solve (rx_rvarf_ / rx_rvar1_<dir> / rx_rvar2_<i>_<j>).  The
#' `rx_..._` naming keeps even constant sensitivities as real output columns, so
#' every piece is read directly from the solve.
#' @noRd
.foceiReadRvar <- function(.d, aug, no) {
  dirs <- aug$dirs; nd <- length(dirs); P2 <- if (is.null(aug$P2r)) aug$P2 else aug$P2r   # rvar2 spans every direction
  aR <- matrix(vapply(dirs, function(q) .d[[paste0("rx_rvar1_", q)]], numeric(no)), no, nd)
  AR <- array(0, c(no, nd, nd))
  for (r in seq_len(nrow(P2))) {
    .ii <- match(P2$i[r], dirs); .jj <- match(P2$j[r], dirs)
    .v2 <- .d[[paste0("rx_rvar2_", P2$i[r], "_", P2$j[r])]]
    AR[, .ii, .jj] <- .v2; AR[, .jj, .ii] <- .v2
  }
  .out <- list(R = .d$rx_rvarf_, aR = aR, AR = AR)
  # sigma pseudo-directions: dR/dsigma (Rsig), d2R/(dsigma ddir) (RsigDir), and
  # d2R/(dsigma dsigma') (Rsig2), keyed by the error-parameter theta indices sigTh.
  .sig <- aug$sigTh
  if (length(.sig) > 0L) {
    .out$Rsig <- matrix(vapply(.sig, function(n) .d[[paste0("rx_rsig_", n, "_")]], numeric(no)), no, length(.sig))
    .out$RsigDir <- array(vapply(.sig, function(n)
      vapply(dirs, function(q) .d[[paste0("rx_rsig1_", n, "_", q)]], numeric(no)), matrix(0, no, nd)),
      c(no, nd, length(.sig)))
    .Rs2 <- array(0, c(no, length(.sig), length(.sig)))
    for (a in seq_along(.sig)) for (b in seq_along(.sig)[seq_along(.sig) >= a]) {
      .v <- .d[[paste0("rx_rsig2_", .sig[a], "_", .sig[b])]]; .Rs2[, a, b] <- .v; .Rs2[, b, a] <- .v
    }
    .out$Rsig2 <- .Rs2
  }
  if (isTRUE(aug$hasTrans))
    .out$trans <- list(yj = .d$rx_tyj_, lambda = .d$rx_tlambda_, low = .d$rx_tlow_, hi = .d$rx_thi_)
  .out
}

#' Re-solve one subject's EBE to the FOCE inner stationary point S_FOCE = sum(q a)
#' + Omega^-1 eta = 0 (q = -eps/R) via Newton on the FOCE inner Hessian
#' Hf = sum(q' a a' + q A) + Omega^-1, starting from the stored eta `eta0`.
#' nlmixr's stored FOCE-combined EBEs do NOT satisfy S_FOCE=0 (an estimation-side
#' inconsistency), so R must be formed at the re-solved eta.  For additive/FOCEI the
#' stored eta is already stationary (|S_FOCE| < `skip`) -> returns `eta0` unchanged
#' (byte no-op).  `NULL` on a solve/Newton failure -> caller falls back to FD.
#' @noRd
.foceiAnalyticFoceEbe <- function(aug, th, eta0, s, times, y, etav, R0, Oi, neta, tol,
                                  maxit = 30L, skip = 1e-3, conv = 1e-9,
                                  foceType = 0L, cens = NULL, limit = NULL) {
  ei <- seq_len(neta)
  # interaction-free FOCE inner gradient/curvature from (f,R0): q0 = -(y-f)/R0 = rho_f,
  # q1 = 1/R0 = rho_ff.  For censored (M2/M3/M4) observations q0/q1 are the EXACT censored
  # rho_f/rho_ff at the frozen R0 (censNormalPartials_) so the re-solved eta* is the censored
  # FOCE stationary point.  foce+ (foceType=1) uses the live conditional R at the trial eta;
  # nonmem freezes R0 at the eta=0 population value passed in.
  .fp <- identical(as.integer(foceType), 1L) || is.null(R0)
  # censored (M2/M3/M4) per-obs CENS + LIMIT (NA/NULL -> uncensored) and the censored-obs index
  .cv <- if (is.null(cens)) integer(length(y)) else as.integer(ifelse(is.na(cens), 0L, cens))
  .lv <- if (is.null(limit)) rep(NA_real_, length(y)) else as.numeric(limit)
  .cw <- which(.cv != 0 | is.finite(.lv))              # censored observations
  .SH <- function(eta) {                               # FOCE S_FOCE and its Jacobian Hf at eta
    E <- .foceiAnalyticSolveFA(aug, c(th, setNames(eta, etav)), s, times, tol = tol)
    if (is.null(E)) return(NULL)
    yt <- .foceiAnalyticTbsY(y, E$trans)               # DV -> rx_pred_ (transformed) scale; no-op if untransformed
    R0e <- if (.fp) E$R else R0
    q0 <- -(yt - E$f) / R0e; q1 <- 1 / R0e
    if (length(.cw)) {                                 # censored: exact rho_f/rho_ff at frozen R0
      .limt <- .foceiAnalyticTbsY(.lv, E$trans)        # transform the censoring bound like the DV
      .cp <- censNormalPartials_(.cv, yt, .limt, E$f, R0e, 2L)
      q0[.cw] <- .cp[.cw, 1]; q1[.cw] <- .cp[.cw, 3]   # cp[,1]=rho_f, cp[,3]=rho_ff
    }
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

#' Batched FOCE/foce+ EBE re-solve: the same interaction-free Newton as
#' [.foceiAnalyticFoceEbe] but over ALL subjects at once via [.foceiAnalyticSolveAll]
#' (one batched solve per Newton iteration instead of per-subject SolveFA).  Bit-identical
#' to the per-subject Newton; avoids the per-subject solve entirely (needed for the shared
#' `dirs` model, which solves batched but not per-subject in the fit's cov-hook context) and
#' is faster.  Returns the nsub x neta eta-hat matrix, or NULL if any subject fails to converge.
#' @noRd
.foceiAnalyticFoceEbeBatch <- function(am, th, ebes, ids, data, obsAll, obsTimes, etav, Oi, neta, tol,
                                       foceType = 0L, E0all = NULL, maxit = 30L, skip = 1e-3, conv = 1e-9) {
  nsub <- nrow(ebes); ei <- seq_len(neta)
  .fp <- identical(as.integer(foceType), 1L) || is.null(E0all)
  Y  <- lapply(obsAll, function(.o) .o$DV)
  CV <- lapply(obsAll, function(.o) if (is.null(.o$CENS)) integer(length(.o$DV)) else as.integer(ifelse(is.na(.o$CENS), 0L, .o$CENS)))
  LV <- lapply(obsAll, function(.o) if (is.null(.o$LIMIT)) rep(NA_real_, length(.o$DV)) else as.numeric(.o$LIMIT))
  .SHi <- function(E, eta_i, i) {                        # S_FOCE + Hf for subject i (censored-aware)
    yt <- .foceiAnalyticTbsY(Y[[i]], E$trans)
    R0e <- if (.fp) E$R else E0all[[i]]$R
    q0 <- -(yt - E$f) / R0e; q1 <- 1 / R0e
    .cw <- which(CV[[i]] != 0 | is.finite(LV[[i]]))
    if (length(.cw)) { .limt <- .foceiAnalyticTbsY(LV[[i]], E$trans)
      .cp <- censNormalPartials_(CV[[i]], yt, .limt, E$f, R0e, 2L); q0[.cw] <- .cp[.cw, 1]; q1[.cw] <- .cp[.cw, 3] }
    S <- as.numeric(Oi %*% eta_i); for (l in ei) S[l] <- S[l] + sum(q0 * E$a[, l])
    Hf <- Oi; for (l in ei) for (m in ei) Hf[l, m] <- Hf[l, m] + sum(q1 * E$a[, l] * E$a[, m] + q0 * E$A[, l, m])
    list(S = S, Hf = Hf)
  }
  eta <- ebes; active <- rep(TRUE, nsub)
  for (it in seq_len(maxit + 1L)) {                      # it=1 evaluates at eta0 (skip test), then Newton steps
    Es <- .foceiAnalyticSolveAll(am, th, eta, ids, data, obsTimes, tol)
    if (is.null(Es)) return(NULL)
    for (i in which(active)) {
      sh <- .SHi(Es[[i]], eta[i, ], i)
      if (max(abs(sh$S)) < (if (it == 1L) skip else conv)) { active[i] <- FALSE; next }
      if (it == maxit + 1L) return(NULL)                 # did not converge -> FD fallback
      step <- tryCatch(solve(sh$Hf, sh$S), error = function(e) NULL); if (is.null(step)) return(NULL)
      eta[i, ] <- eta[i, ] - step
    }
    if (!any(active)) break
  }
  if (any(active)) return(NULL)
  eta
}

#' Compute the full analytic FOCEI covariance (theta + sigma + Omega) for a fitted
#' object, or `NULL` when out of scope / the augmented solve fails.  The cached,
#' env-installing entry point is [foceiCovAnalytic]; this is the raw compute.
#' @param fit a fitted nlmixr2 focei object
#' @return list(cov, se, R, params, method) or `NULL`
#' @noRd
.foceiCovAnalyticCalc <- function(fit) {
  ui <- fit$finalUi
  if (!.hasRxSens())
    return(.foceiAnalyticFallback("an rxode2 without symbolic sensitivities (needs rxExpandSens2_ + symengine)"))
  # FO/FOI is out of scope (the analytic (f,R) path is a FOCEI/FOCE observed information).  The
  # runtime `fo` flag is not persisted to fit$finalUi, so this standalone entry also keys on the
  # persisted estimation method (ui$control$est) -- otherwise an FO fit would be silently assembled
  # as FOCE and mislabelled "analytic".  (The live covType="analytic" hook reads the in-fit ui
  # where `fo` is set, and FO forces covMethod=0, so the production path is already safe.)
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE))) ||
      isTRUE(rxode2::rxGetControl(ui, "est", "") %in% c("fo", "foi")))
    return(.foceiAnalyticFallback("the FO/FOI method"))
  # linCmt() has no symbolic state sensitivities for the augmented model
  if (isTRUE(any(ui$predDf$linCmt)))
    return(.foceiAnalyticFallback("a linCmt() model"))
  interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))                   # 1 FOCEI / 0 FOCE
  # FOCE variance mode (0 "nonmem" frozen R0, 1 "foce+" live R); FOCEI ignores it
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(ui, "foceType", 0L)) else 0L
  if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L)
    return(.foceiAnalyticFallback("adaptive Gaussian quadrature (nAGQ > 1)"))
  # IOV is out of scope; derive the flag from THIS fit, not the process-global
  # .uiIovEnv (which reflects the LAST-preprocessed model).  An occasion (IOV) eta
  # carries a non-"id" `condition` -- the same signal .uiApplyIov keys on.
  .idf0 <- ui$iniDf
  if (any(!is.na(.idf0$condition) & .idf0$condition != "id" & is.na(.idf0$err)))            # IOV
    return(.foceiAnalyticFallback("inter-occasion variability (IOV)"))
  # censored (M2/M3/M4): FOCEI and FOCE with censOption="gauss" are in scope (censored score
  # partials + Gauss-Newton determinant); only the laplace censored determinant uses FD.
  .hasCens <- (!is.null(fit$dataSav$CENS) && any(fit$dataSav$CENS != 0, na.rm = TRUE)) ||
    (!is.null(fit$dataSav$LIMIT) && any(is.finite(fit$dataSav$LIMIT)))
  if (.hasCens && as.integer(rxode2::rxGetControl(ui, "censOption", 0L)) == 1L)
    return(.foceiAnalyticFallback("censored observations with censOption='laplace'"))
  ef <- .foceiAnalyticErrFull(ui)
  if (is.null(ef)) return(NULL)                     # unsupported error model -> errFull already messaged

  ini <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)                    # theta <-> eta pairing
  etaNames <- .map$etaNames
  neta <- length(etaNames)
  if (neta == 0L) return(.foceiAnalyticFallback("a model with no random effects"))
  thetaForEta <- .map$thetaForEta
  # a non-mu-ref (orphan) eta is in scope: it keeps its own ETA_i_ sensitivity direction
  # and an eta-named Omega variance (matching the production hook); no theta maps to it.
  if (any(.iniIsFixed(ini, thetaForEta)))          # fixed structural theta breaks eta indexing
    return(.foceiAnalyticFallback("a fixed mu-referenced structural parameter"))
  keep <- !.iniIsFixed(ini, ef$sgName); ef$sgVar <- ef$sgVar[keep]; ef$sgName <- ef$sgName[keep]  # drop fixed sigma
  Om <- fit$omega
  pairs <- .foceiOmegaPairs(Om, ini)               # free Omega lower-triangle (declared blocks)
  omd <- .omegaVarCovDeriv(Om, pairs)

  # Uniform direction assembly (shared with the production covType="analytic" hook).
  .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta,
                                   sharedEta = unname(.foceiEtaOccurrence(ui) > 1L))
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

  # FOCEI with a general (non-add/prop, multi-endpoint, or estimated-lambda) variance uses the
  # (f,R) cov (sigmas as directions), matching the live covType="analytic" hook; add/prop keeps
  # the fast symbolic assembly.  (IOV already bowed out above.)  An estimated boxCox/yeoJohnson
  # lambda sits in thStruct as a theta-like direction, so name the sigma block with the
  # lambda-excluded .dir$sgName (matching the live hook's fullNm).  Censored FOCEI also routes
  # here (the fast add/prop assembler has no censored partials).
  R <- if (isTRUE(ef$foceiOnly) || .hasCens)
    .foceiAnalyticAssembleRFR(ui, th, ebes, fit$eta$ID, fit$dataSav, Om, ef, neta,
                              length(.dir$dirP), .dir$dirP, omd,
                              dirsCov = .dir$dirsCov, ndirCov = .dir$ndirCov,
                              solveTol = .foceiAnalyticSolveTol(ui),
                              interaction = interaction, foceType = foceType, lamDir = .dir$lamDir)
  else .foceiAnalyticAssembleR(ui, th, ebes, fit$eta$ID, fit$dataSav, Om, ef, neta, nth, nsg, omd,
                               dirs = dirs, dirTh = dirTh, ndir = ndir,
                               solveTol = .foceiAnalyticSolveTol(ui), interaction = interaction,
                               foceType = foceType)
  if (is.null(R)) return(NULL)
  cov <- tryCatch(solve(R), error = function(e) NULL)
  if (is.null(cov)) return(NULL)
  onm <- etaNames                                            # Omega named by the eta (om.eta.cl)
  nm <- c(thStruct, .dir$sgName, .foceiOmegaCovNames(pairs, onm))
  dimnames(R) <- dimnames(cov) <- list(nm, nm)
  list(cov = cov, se = setNames(suppressWarnings(sqrt(diag(cov))), nm),  # NaN flags non-PD
       R = R, params = nm, method = "analytic")
}

#' Full analytic FOCEI covariance (theta + sigma + Omega) for a fitted object, or
#' `NULL` when out of scope / the augmented solve fails.
#'
#' The result is cached on the fit environment (`.covAnalytic`) and its `$cov` is
#' installed as the fit's `$cov` on the first call, so repeated calls -- and
#' `getVarCov()` -- return the stored covariance instead of recomputing the
#' augmented sensitivity solve every time.
#' @param fit a fitted nlmixr2 focei object
#' @return list(cov, se, R, params, method) or `NULL`
#' @noRd
foceiCovAnalytic <- function(fit) {
  .env <- fit
  if (rxode2::rxIs(fit, "nlmixr2FitData")) .env <- fit$env
  if (exists(".covAnalytic", envir = .env, inherits = FALSE)) {
    return(get(".covAnalytic", envir = .env))
  }
  # Match the live covType="analytic" hook (.foceiCalcRanalytic), which wraps the whole assembly
  # in tryCatch and returns NULL on any error -> FD fallback.  A direct foceiCovAnalytic()/
  # getVarCov() call must fall back just as gracefully (e.g. a pure-proportional FOCE fit whose
  # near-zero-prediction branch can hit an NA), never throw.
  .ret <- tryCatch(.foceiCovAnalyticCalc(fit), error = function(e) NULL)
  assign(".covAnalytic", .ret, envir = .env)   # cache (incl. NULL) -- do not recompute
  if (!is.null(.ret) && is.matrix(.ret$cov)) {
    .env$cov <- .ret$cov                        # install so getVarCov()/$cov reuse it
    .env$covMethod <- "analytic"                # report the analytic observed information
  }
  .ret
}
