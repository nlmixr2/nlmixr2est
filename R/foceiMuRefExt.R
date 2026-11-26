#' Extended mu-referencing: eta-sensitivity reuse through the individual-parameter chain rule.
#'
#' Mu-referencing lets a theta reuse its eta's sensitivities when the two enter the model
#' identically (`df/dtheta = df/deta`), and `.foceiAnalyticCovDirMap()` extends this to a
#' covariate coefficient that enters linearly (`df/db = cov*df/deta`).  Both are special
#' cases of one identity: when a theta enters the model only through the individual
#' parameter that carries an eta, on the scale where that eta is additive
#' (`mu(theta, cov) + eta`), then for EVERY model expression E (ODE right-hand side,
#' initial condition, dosing modifier, prediction, residual variance)
#'
#'     dE/dtheta = c * dE/deta,     c = dmu/dtheta,
#'
#' with `c` free of the random effects, the states and time.  Every state sensitivity of
#' the theta is then `c` times the eta's, so the theta needs no sensitivity equations of
#' its own, whatever the (nonlinear) form of mu: maturation (`PMA^hill/(PMA^hill +
#' tm50^hill)`), Emax, natural-scale shifts (`tcl*(1 + b*(WT - 70))`), power exponents
#' outside the exp (`(WT/70)^b`), or a typical value on the natural scale
#' (`tcl*exp(eta.cl)`, c = 1/tcl).
#'
#' Unlike a covariate scale, `c` may depend on thetas, so the 2nd-order columns gain the
#' gradient-times-curvature term `(dc_i/dtheta_j) * f1_eta` (emitted by the caller).
#'
#' Detection is symbolic on the pruned sensitivity environment `.s`:
#'   * candidate: a structural theta with its own direction in `fDirs`, neither reused already
#'     (`covMap`) nor the source another reused direction scales;
#'   * c is the ratio dE/dtheta / dE/deta on an expression where both are nonzero, with the
#'     random effects set to 0 and the states to 1 (an exact simplification when the
#'     identity holds, since c is then free of both);
#'   * c may contain only subject-constant covariates (`constCovs`), fixed thetas, and thetas
#'     with their own direction (so every curvature term (dc_i/dtheta_j) * f1_eta has a column
#'     to land in);
#'   * the identity is then proven on EVERY expression: `expand(dE/dtheta - c*dE/deta)` must be
#'     literally 0.  No numeric fallback: a difference symengine cannot reduce leaves the theta
#'     on its own direction, as does any error -- the full build is always correct.
#'   * models with delay() histories are skipped (the history sensitivities are emitted per
#'     integrated direction by rxode2 and are not part of the checked expressions).
#'
#' Returns a named list `THETA_j_ -> list(etaDir, scale, scaleSE)` in the
#' `.foceiAnalyticCovDirMap()` format plus the symengine scale.  Disabled by
#' `FOCEI_NO_MUREF_EXT=1`.
#' @noRd
.foceiAnalyticMuRefExtDirMap <- function(ui, .s, .st, fDirs, covMap, constCovs) {
  if (nzchar(Sys.getenv("FOCEI_NO_MUREF_EXT"))) {
    return(list())
  }
  # subject-level etas only: an occasion (IOV) eta varies within the subject
  .ini <- ui$iniDf
  .occEta <- .ini$neta1[!is.na(.ini$neta1) & .ini$neta1 == .ini$neta2 &
    !is.na(.ini$condition) & .ini$condition != "id"]
  .etaDirs <- setdiff(grep("^ETA_[0-9]+_$", fDirs, value = TRUE), paste0("ETA_", .occEta, "_"))
  # not already reused, and not the source a reused covariate direction scales (chained reuse)
  .covSrc <- vapply(covMap, function(.m) .m$etaDir, character(1))
  .thDirs <- setdiff(grep("^THETA_[0-9]+_$", fDirs, value = TRUE), c(names(covMap), .covSrc))
  if (length(.etaDirs) == 0L || length(.thDirs) == 0L) {
    return(list())
  }
  # delay() pre-histories carry their own per-direction sensitivities (rxode2 ..pastLines)
  if (!is.null(.s$..pastLines) || length(grep("^rx__pastRhs_", ls(envir = .s, all.names = TRUE)))) {
    return(list())
  }
  .S <- symengine::S
  .get <- function(.v) tryCatch(.S(get(.v, envir = .s)), error = function(e) NULL)
  # every expression the augmented model reads (same sources as .foceiAnalyticAugModelDirs)
  .exprNames <- c(
    paste0("rx__d_dt_", .st, "__"),
    paste0("rx_", .st, "_ini_0__"),
    .foceiAnalyticDoseVars(.s),
    "rx_pred_",
    "rx_r_",
    # transform-of-both-sides parameters, emitted into the augmented model as well
    "rx_yj_", "rx_lambda_", "rx_low_", "rx_hi_"
  )
  .E <- Filter(Negate(is.null), lapply(stats::setNames(.exprNames, .exprNames), .get))
  if (length(.E) == 0L) {
    return(list())
  }
  .isZero <- .rxFoceiCarryIsZero
  .allEta <- unique(unlist(lapply(.E, function(.x) grep("^ETA_[0-9]+_$", .rxFoceiCarryFreeSyms(.x), value = TRUE))))
  # c is free of etas and states when the identity holds, so evaluating the ratio at
  # eta = 0, states = 1 only simplifies it
  .zeroEta <- function(.x) {
    for (.e in .allEta) .x <- symengine::subs(.x, .S(.e), .S(0))
    for (.k in .st) .x <- symengine::subs(.x, .S(.k), .S(1))
    .x
  }
  # symengine keeps symbols that cancel inside a function argument (exp(a + b - (a))); a
  # symbol outside the allowed set is dropped only when c provably does not depend on it
  .dropInert <- function(.x) {
    for (.v in setdiff(.rxFoceiCarryFreeSyms(.x), .okSyms)) {
      if (!.isZero(symengine::expand(symengine::D(.x, .S(.v))))) {
        return(.x) # genuine dependence: the caller's symbol check rejects it
      }
      .x <- symengine::subs(.x, .S(.v), .S(0))
    }
    .x
  }
  .fixedTh <- paste0("THETA_", .ini$ntheta[!is.na(.ini$ntheta) & .ini$fix], "_")
  .okSyms <- c(grep("^THETA_[0-9]+_$", fDirs, value = TRUE), .fixedTh, constCovs)
  # eta derivatives once per eta (shared by every candidate theta)
  .dEt <- lapply(stats::setNames(.etaDirs, .etaDirs), function(.eta) {
    .etS <- .S(.eta)
    lapply(.E, function(.x) symengine::D(.x, .etS))
  })
  .zEt <- lapply(.dEt, function(.d) vapply(.d, .isZero, logical(1)))
  .res <- list()
  for (.th in .thDirs) {
    .thS <- .S(.th)
    .dTh <- lapply(.E, function(.x) symengine::D(.x, .thS))
    .zTh <- vapply(.dTh, .isZero, logical(1))
    if (all(.zTh)) next
    for (.eta in .etaDirs) {
      # an expression where the theta acts but the eta does not breaks the identity
      if (any(!.zTh & .zEt[[.eta]])) next
      .both <- which(!.zTh & !.zEt[[.eta]]) # expressions where both derivatives are nonzero
      # each such expression proposes c; the first proposal proven on all of them is used (an
      # expression that is singular at eta = 0, states = 1 just proposes nothing usable)
      .c <- NULL
      for (.w in .both) {
        .cw <- tryCatch(.dropInert(.zeroEta(.dTh[[.w]] / .dEt[[.eta]][[.w]])), error = function(e) NULL)
        if (is.null(.cw) || .isZero(.cw) || !all(.rxFoceiCarryFreeSyms(.cw) %in% .okSyms)) next
        # proven on every expression the eta moves (where theta does not, c*dE/deta must vanish)
        .proven <- all(vapply(which(!.zEt[[.eta]]), function(.i) {
          .d <- tryCatch(symengine::expand(.dTh[[.i]] - .cw * .dEt[[.eta]][[.i]]), error = function(e) NULL)
          !is.null(.d) && .isZero(.d)
        }, logical(1)))
        if (.proven) {
          .c <- .cw
          break
        }
      }
      if (is.null(.c)) next
      # A scale that is non-finite at some covariate value occurs only where the individual
      # parameter itself is 0 or infinite (dE/deta = 0 there); rxode2's guarded log() keeps
      # such points finite, e.g. c = log(AGE/40) at AGE = 0.
      .res[[.th]] <- list(etaDir = .eta, scale = rxode2::rxFromSE(.c), scaleSE = .c)
      break
    }
  }
  .res
}

#' Gradient-times-curvature term of a 2nd-order column for extended mu-referenced directions.
#'
#' For a pair (i, j) in which direction i reuses eta e_i with a theta-dependent scale c_i,
#'   d2f/(di dj) = c_i * c_j * f2[e_i, e_j] + (dc_i/dtheta_j) * f1[e_i],
#' (c_j = 1 and e_j = j when j is not scaled).  The product part is emitted by the caller;
#' this returns the second term as rxode2 text ("" when it vanishes).  Only one of the two
#' orderings is used: if both are scaled, symmetry of d2mu/(dtheta_i dtheta_j) makes the
#' two candidate terms equal.
#' @noRd
.foceiMuRefExtCurvTerm <- function(map, i, j, pre1) {
  .one <- function(.a, .b) {
    .m <- map[[.a]]
    if (is.null(.m) || is.null(.m$scaleSE) || !grepl("^THETA_[0-9]+_$", .b)) {
      return("")
    }
    .d <- symengine::D(.m$scaleSE, symengine::S(.b))
    if (.rxFoceiCarryIsZero(.d)) {
      return("")
    }
    paste0("+(", rxode2::rxFromSE(.d), ")*", pre1, .m$etaDir)
  }
  .t <- .one(i, j)
  if (nzchar(.t)) .t else .one(j, i)
}
