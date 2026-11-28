#' Extended mu-referencing: thetas whose sensitivities are a subject-constant multiple of an eta's.
#'
#' A theta acting only through an eta-carrying parameter satisfies dE/dtheta = c * dE/deta for
#' every model expression E, with c = dmu/dtheta free of etas, states and time; its columns are
#' then emitted as c times the eta's (the identity is proven symbolically, else no reuse).
#' @param ui rxode2 ui
#' @param .s pruned sensitivity environment (`ui$loadPruneSens`)
#' @param .st ODE state names
#' @param fDirs model directions of the augmented model
#' @param covMap covariate-coefficient reuse map (`.foceiAnalyticCovDirMap()`)
#' @param constCovs subject-constant covariates
#' @return named list `THETA_j_ -> list(etaDir, scale, scaleSE)`; empty when disabled
#'   (any non-empty `FOCEI_NO_MUREF_EXT`) or nothing qualifies
#' @noRd
.foceiAnalyticMuRefExtDirMap <- function(ui, .s, .st, fDirs, covMap, constCovs) {
  # occasion etas (rewritten as subject-level by the IOV expansion) vary within the subject
  if (nzchar(Sys.getenv("FOCEI_NO_MUREF_EXT")) || .uiHasIov(ui)) {
    return(list())
  }
  .ini <- ui$iniDf
  # error-model thetas would need an rsig1 curvature term
  .errTh <- paste0("THETA_", .ini$ntheta[!is.na(.ini$ntheta) & !is.na(.ini$err)], "_")
  .etaDirs <- grep("^ETA_[0-9]+_$", fDirs, value = TRUE)
  .thAll <- setdiff(grep("^THETA_[0-9]+_$", fDirs, value = TRUE), .errTh)
  # not already reused, and not the source theta of an eta-less covariate entry: that entry's
  # column reads the source's, which a reused source would only emit later
  .covSrc <- vapply(covMap, function(.m) .m$etaDir, character(1))
  .thDirs <- setdiff(.thAll, c(names(covMap), .covSrc))
  if (length(.etaDirs) == 0L || length(.thDirs) == 0L) {
    return(list())
  }
  .vars <- ls(envir = .s, all.names = TRUE)
  .odeNames <- paste0("rx__d_dt_", .st, "__")
  # matExp()/indLin() ODEs are only materialized later (.rxJacobian); delay() histories carry
  # their own per-direction sensitivities; mtime() right-hand sides are not part of the proof
  if (!all(.odeNames %in% .vars) || any(startsWith(.vars, "rx__pastRhs_")) || nzchar(.mtimeLinesStr(.s))) {
    return(list())
  }
  .S <- symengine::S
  .get <- function(.v) tryCatch(.S(get(.v, envir = .s)), error = function(e) NULL)
  .exprNames <- c(
    .odeNames,
    paste0("rx_", .st, "_ini_0__"),
    .foceiAnalyticDoseVars(.s),
    "rx_pred_",
    "rx_r_",
    "rx_yj_",
    "rx_lambda_",
    "rx_low_",
    "rx_hi_"
  )
  .exprNames <- intersect(.exprNames, .vars)
  .E <- lapply(stats::setNames(.exprNames, .exprNames), .get)
  # an expression symengine cannot represent must not silently drop out of the proof
  if (
    length(.E) == 0L ||
      any(vapply(.E, is.null, logical(1))) ||
      any(vapply(.E, function(.x) grepl("delay(", rxode2::rxFromSE(.x), fixed = TRUE), logical(1)))
  ) {
    return(list())
  }
  .isZero <- .rxFoceiCarryIsZero
  .allEta <- unique(unlist(lapply(.E, function(.x) grep("^ETA_[0-9]+_$", .rxFoceiCarryFreeSyms(.x), value = TRUE))))
  # c is free of etas and states when the identity holds, so these substitutions only simplify it
  .zeroEta <- function(.x) {
    for (.e in .allEta) {
      .x <- symengine::subs(.x, .S(.e), .S(0))
    }
    for (.k in .st) {
      .x <- symengine::subs(.x, .S(.k), .S(1))
    }
    .x
  }
  .fixedTh <- paste0("THETA_", .ini$ntheta[!is.na(.ini$ntheta) & .ini$fix], "_")
  .okSyms <- c(.thAll, .fixedTh, constCovs)
  # symengine keeps symbols that cancel inside a function argument; drop one only when dc/dsym == 0
  .dropInert <- function(.x) {
    for (.v in setdiff(.rxFoceiCarryFreeSyms(.x), .okSyms)) {
      if (!.isZero(symengine::expand(symengine::D(.x, .S(.v))))) {
        return(.x)
      }
      .x <- symengine::subs(.x, .S(.v), .S(0))
    }
    .x
  }
  .dEt <- lapply(stats::setNames(.etaDirs, .etaDirs), function(.eta) {
    .etS <- .S(.eta)
    lapply(.E, function(.x) symengine::D(.x, .etS))
  })
  .zEt <- lapply(.dEt, function(.d) vapply(.d, .isZero, logical(1)))
  # cheap rejection before expand(): a clearly nonzero value at one point disproves the identity
  .clearlyNonzero <- function(.d, .ref) {
    tryCatch(
      {
        .sy <- union(.rxFoceiCarryFreeSyms(.d), .rxFoceiCarryFreeSyms(.ref))
        for (.k in seq_along(.sy)) {
          .val <- .S(0.31 + (.k * 0.6180339887) %% 1)
          .d <- symengine::subs(.d, .S(.sy[.k]), .val)
          .ref <- symengine::subs(.ref, .S(.sy[.k]), .val)
        }
        .vd <- as.double(symengine::evalf(.d))
        .vr <- as.double(symengine::evalf(.ref))
        is.finite(.vd) && is.finite(.vr) && abs(.vd) > 1e-6 * (abs(.vr) + 1)
      },
      error = function(e) FALSE
    )
  }
  .res <- list()
  for (.th in .thDirs) {
    .thS <- .S(.th)
    .dTh <- lapply(.E, function(.x) symengine::D(.x, .thS))
    .zTh <- vapply(.dTh, .isZero, logical(1))
    if (all(.zTh)) {
      next
    }
    for (.eta in .etaDirs) {
      if (any(!.zTh & .zEt[[.eta]])) {
        next
      } # theta acts where the eta does not
      .c <- NULL
      .tried <- character(0)
      for (.w in which(!.zTh & !.zEt[[.eta]])) {
        .cw <- tryCatch(.dropInert(.zeroEta(.dTh[[.w]] / .dEt[[.eta]][[.w]])), error = function(e) NULL)
        # a scale whose symbols cannot be listed is refused, not assumed constant
        .fs <- if (is.null(.cw)) {
          NULL
        } else {
          tryCatch(vapply(symengine::free_symbols(.cw), as.character, ""), error = function(e) NULL)
        }
        if (is.null(.fs) || .isZero(.cw) || !all(.fs %in% .okSyms)) {
          next
        }
        # the same candidate from another expression would repeat the same proof
        .key <- as.character(.cw)
        if (.key %in% .tried) {
          next
        }
        .tried <- c(.tried, .key)
        # proven on every expression the eta moves
        .idx <- which(!.zEt[[.eta]])
        .diff <- lapply(.idx, function(.i) .dTh[[.i]] - .cw * .dEt[[.eta]][[.i]])
        if (any(mapply(.clearlyNonzero, .diff, lapply(.idx, function(.i) .dEt[[.eta]][[.i]])))) {
          next
        }
        .proven <- all(vapply(
          .diff,
          function(.d) {
            .d <- tryCatch(symengine::expand(.d), error = function(e) NULL)
            !is.null(.d) && .isZero(.d)
          },
          logical(1)
        ))
        if (.proven) {
          .c <- .cw
          break
        }
      }
      if (is.null(.c)) {
        next
      }
      .res[[.th]] <- list(etaDir = .eta, scale = rxode2::rxFromSE(.c), scaleSE = .c)
      break
    }
  }
  .res
}

#' Curvature term (dc_i/dtheta_j) * f1[e_i] of a 2nd-order pair column with a theta-dependent scale.
#'
#' Either ordering gives the same term (symmetry of the second derivative of mu), so the first
#' nonzero one is used.
#' @param map reuse map (`covMap` merged with the extended entries)
#' @param i,j the pair's direction names
#' @param pre1 1st-order column prefix (`rx_f1_` or `rx_rvar1_`)
#' @return rxode2 text to append, or ""
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
