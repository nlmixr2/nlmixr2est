# Analytic FOCEI/FOCE observed-information covariance R-matrix for nlmixr2est.
#
# Replaces the finite-difference Hessian (`foceiCalcR`) with the exact analytic
# observed information -- the exact-gradient construction continued by one
# differentiation.  The population Hessian splits into a data term (reusing the
# 2nd-order sensitivities already needed for the gradient) and a log-determinant
# term (which carries the 3rd-order sensitivities), assembled from 1st/2nd/3rd
# -order model sensitivities (rxExpandSens_/2_/3_ in rxode2).  Validated against
# ferx (FeRx-NLME/ferx-core#450) and finite differences of the objective.
#
# ROBUSTNESS (per design + ferx all-or-nothing fallback): the augmented
# sensitivity ODE is large and more likely to fail solving than the plain
# model.  Every per-subject solve is guarded; ANY failure (solve error,
# non-finite sensitivity, out-of-scope model, or an rxode2 that lacks
# `rxExpandSens3_`) makes the whole analytic path return `NULL`, and the caller
# (foceiCov) drops to the next tier (fd2, then the builtin objective FD).  The
# analytic path is therefore never partially applied.

#' Install the full natural-scale analytic covariance as the fit's `$cov`.
#'
#' When `foceiCalcR` -> `.foceiCalcRanalytic` ran (analytic exact3 or fd2 tier),
#' the full theta + sigma + Omega natural-scale covariance is stashed as
#' `.analyticCov`.  The native cov path only fills the theta block of `$cov`
#' (Omega/residual are skipCov'd), so here we replace `$cov` with the full
#' matrix.  Users extract Omega/residual SEs themselves via `sqrt(diag(fit$cov))`.
#' On the finite-difference fallback (no `.analyticCov`) `$cov` stays theta-only.
#' @param .ret focei fit environment
#' @return Nothing; mutates `.ret$cov`
#' @noRd
.foceiInstallAnalyticCov <- function(.ret) {
  if (!exists(".analyticCov", envir = .ret, inherits = FALSE)) return(invisible())
  .cov <- get(".analyticCov", envir = .ret)
  if (!is.matrix(.cov) || !all(is.finite(.cov))) return(invisible())
  .ret$cov <- .cov                       # analytic-tier cov already carries dimnames
  invisible()
}

#' Native-path hook: build the analytic observed-information R-matrix for the
#' running focei fit, called from C++ `foceiCalcR` for an R-matrix covariance.
#'
#' Sources the converged UI / data / Omega / EBEs / estimates from the live focei
#' environment `e`, builds the augmented sensitivity model, and assembles the full
#' natural-scale observed-information R (theta + sigma + Omega) per subject via
#' [.foceiAnalyticSubjectR] (Omega block from [.omegaVarCovDeriv] /
#' `rxOmegaVarCovDeriv`), summed over subjects.  Tries the exact 3rd-order tier,
#' then the fd2 tier; stashes the full natural cov (installed as `fit$cov` by
#' [.foceiInstallAnalyticCov]) and returns the theta `R.0` for the native theta-SE
#' path, or `NULL` to fall back to the finite-difference Hessian.
#' @noRd
.foceiCalcRanalytic <- function(e) {
  tryCatch({
    if (!.hasRxExpandSens3()) return(NULL)
    ui <- get("ui", e)
    ef <- .foceiAnalyticErrFull(ui)
    if (is.null(ef)) return(NULL)

    ini <- ui$iniDf; muRef <- ui$muRefDataFrame
    etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
    etaRows <- etaRows[order(etaRows$neta1), , drop = FALSE]
    etaNames <- etaRows$name; neta <- length(etaNames)
    if (neta == 0L) return(NULL)
    thetaForEta <- muRef$theta[match(etaNames, muRef$eta)]
    if (anyNA(thetaForEta)) return(NULL)             # non-mu-ref eta -> out of scope

    # converged estimates (e$theta, NOT ui$iniDf$est which holds the initials)
    thNames <- get("thetaNames", e)
    thVals  <- get("theta", e)$theta
    names(thVals) <- thNames
    # cov parameters = non-skipped thetas (skipCov drops residual/fixed/IOV)
    skip <- get("skipCov", e)[seq_along(thNames)]
    covParams <- thNames[!skip]
    if (length(setdiff(covParams, thetaForEta)) > 0L) return(NULL)  # only mu-ref structural thetas
    idx <- match(covParams, thetaForEta)             # eta-order -> cov-param order
    if (anyNA(idx)) return(NULL)

    # converged residual sigma into the error-model evaluator (ef reads initials)
    .valc <- setNames(as.numeric(thVals[ef$sgName]), ef$sgVar)
    ef$ev <- local({ v <- .valc; function(expr, f, y) eval(expr, c(list(f = f, y = y), as.list(v))) })

    # scope: fixed structural thetas break the eta<->theta indexing
    fx <- function(nm) { if (is.null(ini$fix)) return(rep(FALSE, length(nm)))
      i <- match(nm, ini$name); !is.na(i) & !is.na(ini$fix[i]) & ini$fix[i] }
    if (any(fx(thetaForEta))) return(NULL)
    keep <- !fx(ef$sgName); ef$sgVar <- ef$sgVar[keep]; ef$sgName <- ef$sgName[keep]  # drop fixed sigma

    Om <- get("omega", e)
    blocks <- .omegaBlocks(Om)                        # free Omega lower-triangle (diag or block)
    pairs <- do.call(rbind, lapply(blocks, function(b)
      do.call(rbind, lapply(seq_along(b), function(a) cbind(b[a], b[seq_len(a)])))))
    pairs <- pairs[!.omegaFixed(ini, pairs), , drop = FALSE]
    omd <- .omegaVarCovDeriv(Om, pairs)               # non-Cholesky Omega derivs (rxOmegaVarCovDeriv)
    nth <- length(thetaForEta); nsg <- length(ef$sgVar)
    omNm <- apply(pairs, 1, function(p) if (p[1] == p[2]) paste0("om.", thetaForEta[p[1]])
                  else paste0("cov.", thetaForEta[p[1]], ".", thetaForEta[p[2]]))
    fullNm <- c(thetaForEta, ef$sgName, omNm)         # full natural-scale parameter order

    th <- setNames(as.numeric(thVals[thNames]), paste0("THETA_", seq_along(thNames), "_"))
    etav <- paste0("ETA_", seq_len(neta), "_")
    etaObf <- get("etaObf", e)
    ebes <- as.matrix(etaObf[, paste0("ETA[", seq_len(neta), "]"), drop = FALSE])
    ids  <- etaObf$ID
    data <- get("dataSav", e)

    # Assemble the FULL natural-scale observed-information R (theta + sigma + Omega)
    # per subject via .foceiAnalyticSubjectR (which uses `omd` for the Omega block).
    .assemble <- function(sens) {
      am <- if (sens == "fd2") .foceiAnalyticSens2Model(ui) else .foceiAnalyticAugModel(ui)
      if (is.null(am) || am$neta != neta) return(NULL)
      np <- nth + nsg + omd$nom
      R <- matrix(0, np, np)
      for (i in seq_along(ids)) {
        s <- data[data$ID == ids[i], , drop = FALSE]
        obs <- s[s$EVID == 0, , drop = FALSE]
        E <- if (sens == "fd2")
          .foceiAnalyticSolveSubjectFD3(am, c(th, setNames(ebes[i, ], etav)), s, obs$TIME, etav, neta, am$P2, ebes[i, ])
        else
          .foceiAnalyticSolveSubject(am$augMod, c(th, setNames(ebes[i, ], etav)), s, obs$TIME, etav, neta, am$P2, am$P3)
        if (is.null(E)) return(NULL)
        E$y <- obs$DV
        Ri <- tryCatch(.foceiAnalyticSubjectR(E, ebes[i, ], Om, ef, neta, nth, nsg, ef$sgVar, omd),
                       error = function(x) NULL)
        if (is.null(Ri) || !all(is.finite(Ri))) return(NULL)
        R <- R + Ri
      }
      R
    }

    Rfull <- .assemble("exact3")
    if (is.null(Rfull)) Rfull <- .assemble("fd2")
    if (is.null(Rfull)) return(NULL)
    dimnames(Rfull) <- list(fullNm, fullNm)

    # full natural-scale covariance (theta + sigma + Omega) -> becomes fit$cov
    # (people extract the Omega/residual SEs from it; theta SEs flow through the
    # native path below).  The native FD fallback covers theta only.
    covFull <- tryCatch(solve(Rfull), error = function(x) NULL)
    if (is.null(covFull)) return(NULL)
    dimnames(covFull) <- list(fullNm, fullNm)
    assign(".analyticCov",    covFull, envir = e)
    assign(".analyticSE",     setNames(suppressWarnings(sqrt(diag(covFull))), fullNm), envir = e)
    assign(".analyticParams", fullNm,  envir = e)

    # R.0 for the native theta path: precision of the theta block of the FULL
    # inverse, so the native theta SEs are consistent with fit$cov (full matrix).
    .ct <- covFull[seq_len(nth), seq_len(nth), drop = FALSE][idx, idx, drop = FALSE]
    solve(.ct)
  }, error = function(x) NULL)
}

#' 3rd-order sensitivities available? Always TRUE now — rxode2 C++ version if
#' present, otherwise the built-in pure-R generator. Kept for the scope gate.
#' @return logical
#' @noRd
.hasRxExpandSens3 <- function() {
  # requires only rxExpandSens2_ + symengine from rxode2 (>= 5.0.2)
  exists("rxExpandSens2_", envir = asNamespace("rxode2"), inherits = FALSE) &&
    requireNamespace("symengine", quietly = TRUE)
}

#' Non-Cholesky Omega derivatives for the analytic Omega block: Omega^{-1} and
#' log|Omega| first/second derivatives w.r.t. the free variance-covariance
#' elements `pairs` (each row `c(a, b)`, `a >= b`).  Uses rxode2's
#' `rxOmegaVarCovDeriv` when present (nlmixr2/rxode2#1092), else the identical
#' E-basis matrix calculus inline, so block Omega runs on released rxode2.
#' Returns `dOi` (list of dOmega^{-1}/dw), `d2Oi` (list-of-lists), `dLD`, `d2LD`.
#' @noRd
.omegaVarCovDeriv <- function(Om, pairs) {
  neta <- nrow(Om); nom <- nrow(pairs)
  if (exists("rxOmegaVarCovDeriv", envir = asNamespace("rxode2"), inherits = FALSE)) {
    d <- rxode2::rxOmegaVarCovDeriv(Om, order = 2L)
    key <- function(P) paste(P[, 1], P[, 2], sep = "-")
    idx <- match(key(pairs), key(d$elements))
    return(list(nom = nom, dOi = d$dOmegaInv[idx],
                d2Oi = lapply(idx, function(a) d$d2OmegaInv[[a]][idx]),
                d2LD = d$d2LogDet[idx, idx, drop = FALSE]))
  }
  Oi <- solve(Om)
  Eb <- lapply(seq_len(nom), function(m) { E <- matrix(0, neta, neta)
    a <- pairs[m, 1]; b <- pairs[m, 2]; E[a, b] <- 1; E[b, a] <- 1; E })
  OiE <- lapply(Eb, function(E) Oi %*% E)
  list(nom = nom,
       dOi = lapply(OiE, function(M) -(M %*% Oi)),
       d2Oi = lapply(seq_len(nom), function(a) lapply(seq_len(nom), function(b)
         Oi %*% (Eb[[a]] %*% OiE[[b]] + Eb[[b]] %*% OiE[[a]]) %*% Oi)),
       d2LD = outer(seq_len(nom), seq_len(nom),
                    Vectorize(function(a, b) -sum(diag(OiE[[a]] %*% OiE[[b]])))))
}

#' Full error machinery for the theta+sigma+Omega R-matrix: the f-derivatives of
#' rho and p (for the eta/theta blocks and H~), AND the residual-parameter (sigma)
#' partials needed for the sigma block.  Built symbolically in (f, y, sa, sp)
#' where sa = additive SD, sp = proportional SD.  Returns `NULL` for out-of-scope
#' error models (combined1 / lnorm -> FD fallback).
#' @noRd
.foceiAnalyticErrFull <- function(ui) {
  ini <- ui$iniDf; er <- ini[!is.na(ini$err), , drop = FALSE]
  if (any(er$err %in% c("lnorm", "propT", "propF"))) return(NULL)
  addPr <- tryCatch(rxode2::rxGetControl(ui, "addProp", "combined2"), error = function(e) "combined2")
  if (identical(addPr, "combined1")) return(NULL)
  addN <- er$name[er$err == "add"]; propN <- er$name[er$err == "prop"]
  hasA <- length(addN) == 1L; hasP <- length(propN) == 1L
  if (!hasA && !hasP) return(NULL)
  Rstr <- if (hasA && hasP) "sa^2+sp^2*f^2" else if (hasP) "sp^2*f^2" else "sa^2"
  Rq <- parse(text = Rstr)[[1]]
  rhoE <- bquote(0.5 * ((y - f)^2 / .(Rq) + log(.(Rq))))
  pE <- bquote(1 / .(Rq) + 0.5 * (.(D(Rq, "f")) / .(Rq))^2)
  DD <- function(e, ...) { for (v in c(...)) e <- D(e, v); e }
  sgVar <- c(if (hasA) "sa", if (hasP) "sp"); sgName <- c(if (hasA) addN, if (hasP) propN)
  val <- setNames(c(if (hasA) er$est[er$name == addN], if (hasP) er$est[er$name == propN]), sgVar)
  sc <- list(r1 = DD(rhoE, "f"), r2 = DD(rhoE, "f", "f"), r3 = DD(rhoE, "f", "f", "f"),
             p = pE, p1 = DD(pE, "f"), p2 = DD(pE, "f", "f"))
  per <- list(); for (s in sgVar) per[[s]] <- list(rf = DD(rhoE, "f", s), rff = DD(rhoE, "f", "f", s),
    ps = DD(pE, s), pf = DD(pE, "f", s))
  pair <- list(); for (i in seq_along(sgVar)) for (j in i:length(sgVar)) { a <- sgVar[i]; b <- sgVar[j]
    pair[[paste0(a, b)]] <- list(rss = DD(rhoE, a, b), rfss = DD(rhoE, "f", a, b), pss = DD(pE, a, b)) }
  list(sgVar = sgVar, sgName = sgName, sc = sc, per = per, pair = pair,
       ev = function(e, f, y) eval(e, c(list(f = f, y = y), as.list(val))))
}

#' Build the augmented rxode2 model carrying df/deta to 3rd order.
#'
#' Reuses the focei symengine pipeline (`ui$foceiEtaS` gives the 1st-order eta
#' state sensitivities + `rx_pred_`); adds 2nd/3rd-order state sensitivities via
#' `rxExpandSens2_`/`rxExpandSens3_` and the prediction chain f1/f2/f3 by total
#' differentiation of `rx_pred_`.  `ETA[k]`/`THETA[k]` are rewritten to plain
#' params so the model can be solved at arbitrary eta/theta.
#'
#' @param ui rxode2 UI object
#' @return list(augMod, etav, neta, st, P2, P3) or `NULL` on any failure
#' @noRd
.foceiAnalyticAugModel <- function(ui) {
  tryCatch({
    s <- ui$foceiEtaS
    st   <- rxode2::rxStateOde(s)
    neta <- s$..maxEta
    etav <- paste0("ETA_", seq_len(neta), "_")
    s1 <- s$..sens
    s2 <- rxode2::.rxSens(s, etav, etav)
    # 3rd-order sensitivities via rxode2's .rxSens (which dispatches to the C++
    # rxExpandSens3_), mirroring the 2nd-order call above -- one convention for
    # every sensitivity order.
    s3 <- rxode2::.rxSens(s, etav, etav, etav)

    .Dn  <- function(e, v) symengine::D(e, symengine::S(v))
    .sn1 <- function(j, ...) symengine::S(paste0("rx__sens_", j, "_BY_",
                                                 paste(c(...), collapse = "_BY_"), "__"))
    .toRx <- function(.l) rxode2::rxFromSE(.l)
    .g1 <- function(ex, p) { e <- .Dn(ex, p); for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p); e }
    .g2 <- function(ex, p, q) { gq <- .g1(ex, q); e <- .Dn(gq, p)
      for (k in st) e <- e + .Dn(gq, k) * .sn1(k, p); for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p, q); e }
    .g3 <- function(ex, p, q, r) { gqr <- .g2(ex, q, r); e <- .Dn(gqr, p)
      for (k in st) e <- e + .Dn(gqr, k) * .sn1(k, p)
      for (aa in unique(c(q, r))) for (m in st) e <- e + symengine::D(gqr, .sn1(m, aa)) * .sn1(m, p, aa)
      for (m in st) e <- e + symengine::D(gqr, .sn1(m, q, r)) * .sn1(m, p, q, r); e }

    pred <- get("rx_pred_", s)
    P2 <- expand.grid(i = etav, j = etav, stringsAsFactors = FALSE)
    P3 <- expand.grid(i = etav, j = etav, k = etav, stringsAsFactors = FALSE)
    fL1 <- vapply(etav, function(p) { .l <- .g1(pred, p); paste0("f1_", p, "=", .toRx(.l)) }, character(1))
    fL2 <- vapply(seq_len(nrow(P2)), function(r) { .l <- .g2(pred, P2$i[r], P2$j[r])
      paste0("f2_", P2$i[r], "_", P2$j[r], "=", .toRx(.l)) }, character(1))
    fL3 <- vapply(seq_len(nrow(P3)), function(r) { .l <- .g3(pred, P3$i[r], P3$j[r], P3$k[r])
      paste0("f3_", P3$i[r], "_", P3$j[r], "_", P3$k[r], "=", .toRx(.l)) }, character(1))
    baseOde <- vapply(st, function(x) {
      .l <- get(paste0("rx__d_dt_", x, "__"), s); paste0("d/dt(", x, ")=", .toRx(.l)) }, character(1))
    modTxt <- paste(c(baseOde, s1, s2, s3, paste0("predf=", .toRx(pred)), fL1, fL2, fL3), collapse = "\n")
    modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", modTxt)
    modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", modTxt)
    list(augMod = rxode2::rxode2(modTxt), etav = etav, neta = neta, st = st, P2 = P2, P3 = P3)
  }, error = function(e) NULL)
}

#' Solve the augmented sensitivity model for one subject, guarding failures.
#'
#' @return list with per-observation f and the eta-sensitivity arrays
#'   (`a`, `A`, `Ath`), or `NULL` on ANY solve failure / non-finite result.
#' @noRd
.foceiAnalyticSolveSubject <- function(augMod, params, ev, times, etav, neta, P2, P3) {
  .d <- tryCatch(
    as.data.frame(rxode2::rxSolve(augMod, params = params, ev,
                                  returnType = "data.frame", atol = 1e-10, rtol = 1e-10)),
    error = function(e) NULL, warning = function(w) NULL)
  if (is.null(.d)) return(NULL)
  .d <- .d[.d$time %in% times, , drop = FALSE]
  if (nrow(.d) == 0L) return(NULL)
  .cols <- c("predf", paste0("f1_", etav))
  if (!all(.cols %in% names(.d))) return(NULL)
  .a <- sapply(etav, function(p) .d[[paste0("f1_", p)]])
  .A <- array(0, c(nrow(.d), neta, neta))
  .At <- array(0, c(nrow(.d), neta, neta, neta))
  for (r in seq_len(nrow(P2)))
    .A[, match(P2$i[r], etav), match(P2$j[r], etav)] <- .d[[paste0("f2_", P2$i[r], "_", P2$j[r])]]
  for (r in seq_len(nrow(P3)))
    .At[, match(P3$i[r], etav), match(P3$j[r], etav), match(P3$k[r], etav)] <-
      .d[[paste0("f3_", P3$i[r], "_", P3$j[r], "_", P3$k[r])]]
  .out <- list(f = .d$predf, a = matrix(.a, nrow(.d), neta), A = .A, Ath = .At)
  # guard: any non-finite sensitivity invalidates the whole analytic path
  if (!all(is.finite(.out$f)) || !all(is.finite(.out$a)) ||
      !all(is.finite(.out$A)) || !all(is.finite(.out$Ath))) return(NULL)
  .out
}

#' Augmented 3rd-order sensitivity model over an ARBITRARY direction set `dirs`
#' (a mix of `ETA_i_` and `THETA_j_`).  Generalizes [.foceiAnalyticAugModel] (which
#' is the special case dirs = etas): every structural theta differentiates in one
#' direction -- a mu-ref theta reuses its eta's direction (free), a non-mu-ref
#' structural theta (covariate coefficient OR eta-less) gets its own `THETA_j_`
#' direction carrying its TRUE sensitivity.  Uses `rxExpandSens_/2_/3_` over `dirs`,
#' so no covariate is "detected" or scaled -- correctness is by construction.
#' @return list(augMod, dirs, ndir, st, P2, P3) or `NULL`
#' @noRd
.foceiAnalyticAugModelDirs <- function(ui, dirs) {
  tryCatch({
    s <- ui$loadPruneSens
    st <- rxode2::rxStateOde(s)
    rxode2::.rxJacobian(s, c(st, dirs))
    mat <- function(grd) {
      model <- s   # the rxExpandSens "line" expressions reference `model`
      lapply(c(grd$ddtS, grd$ddS2), function(x) assign(x, symengine::Symbol(x), envir = s))
      apply(grd, 1, function(x) { .l <- eval(parse(text = x[["line"]])); paste0(x[["ddt"]], "=", rxode2::rxFromSE(.l)) })
    }
    s1 <- mat(rxode2::rxExpandSens_(st, dirs))
    s2 <- mat(rxode2::rxExpandSens2_(st, dirs, dirs))
    s3 <- mat(rxode2::rxExpandSens3_(st, dirs, dirs, dirs))
    pred <- get("rx_pred_", s)
    .Dn  <- function(e, v) symengine::D(e, symengine::S(v))
    .sn1 <- function(j, ...) symengine::S(paste0("rx__sens_", j, "_BY_", paste(c(...), collapse = "_BY_"), "__"))
    .toRx <- function(.l) rxode2::rxFromSE(.l)
    .g1 <- function(ex, p) { e <- .Dn(ex, p); for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p); e }
    .g2 <- function(ex, p, q) { gq <- .g1(ex, q); e <- .Dn(gq, p)
      for (k in st) e <- e + .Dn(gq, k) * .sn1(k, p); for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p, q); e }
    .g3 <- function(ex, p, q, r) { gqr <- .g2(ex, q, r); e <- .Dn(gqr, p)
      for (k in st) e <- e + .Dn(gqr, k) * .sn1(k, p)
      for (aa in unique(c(q, r))) for (m in st) e <- e + symengine::D(gqr, .sn1(m, aa)) * .sn1(m, p, aa)
      for (m in st) e <- e + symengine::D(gqr, .sn1(m, q, r)) * .sn1(m, p, q, r); e }
    P2 <- expand.grid(i = dirs, j = dirs, stringsAsFactors = FALSE)
    P3 <- expand.grid(i = dirs, j = dirs, k = dirs, stringsAsFactors = FALSE)
    fL1 <- vapply(dirs, function(p) paste0("f1_", p, "=", .toRx(.g1(pred, p))), character(1))
    fL2 <- vapply(seq_len(nrow(P2)), function(r) paste0("f2_", P2$i[r], "_", P2$j[r], "=", .toRx(.g2(pred, P2$i[r], P2$j[r]))), character(1))
    fL3 <- vapply(seq_len(nrow(P3)), function(r) paste0("f3_", P3$i[r], "_", P3$j[r], "_", P3$k[r], "=", .toRx(.g3(pred, P3$i[r], P3$j[r], P3$k[r]))), character(1))
    baseOde <- vapply(st, function(x) paste0("d/dt(", x, ")=", .toRx(get(paste0("rx__d_dt_", x, "__"), s))), character(1))
    modTxt <- paste(c(baseOde, s1, s2, s3, paste0("predf=", .toRx(pred)), fL1, fL2, fL3), collapse = "\n")
    modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", modTxt)
    modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", modTxt)
    list(augMod = rxode2::rxode2(modTxt), dirs = dirs, ndir = length(dirs), st = st, P2 = P2, P3 = P3)
  }, error = function(e) NULL)
}

#' Solve the direction-set augmented model for one subject; returns per-obs `f` and
#' the direction-indexed sensitivity arrays `a` (no x ndir), `A` (no x ndir x ndir),
#' `Ath` (no x ndir^3).  `NULL` on any solve failure / non-finite result.
#' @noRd
.foceiAnalyticSolveDir <- function(aug, params, ev, times) {
  .d <- tryCatch(as.data.frame(rxode2::rxSolve(aug$augMod, params = params, ev,
                                               returnType = "data.frame", atol = 1e-10, rtol = 1e-10)),
                 error = function(e) NULL, warning = function(w) NULL)
  if (is.null(.d)) return(NULL)
  .d <- .d[.d$time %in% times, , drop = FALSE]
  if (nrow(.d) == 0L) return(NULL)
  dirs <- aug$dirs; nd <- length(dirs); no <- nrow(.d)
  if (!all(c("predf", paste0("f1_", dirs)) %in% names(.d))) return(NULL)
  .a <- matrix(vapply(dirs, function(p) .d[[paste0("f1_", p)]], numeric(no)), no, nd)
  .A <- array(0, c(no, nd, nd))
  for (r in seq_len(nrow(aug$P2)))
    .A[, match(aug$P2$i[r], dirs), match(aug$P2$j[r], dirs)] <- .d[[paste0("f2_", aug$P2$i[r], "_", aug$P2$j[r])]]
  .At <- array(0, c(no, nd, nd, nd))
  for (r in seq_len(nrow(aug$P3)))
    .At[, match(aug$P3$i[r], dirs), match(aug$P3$j[r], dirs), match(aug$P3$k[r], dirs)] <-
      .d[[paste0("f3_", aug$P3$i[r], "_", aug$P3$j[r], "_", aug$P3$k[r])]]
  .out <- list(f = .d$predf, a = .a, A = .A, Ath = .At)
  if (!all(is.finite(.out$f)) || !all(is.finite(.out$a)) ||
      !all(is.finite(.out$A)) || !all(is.finite(.out$Ath))) return(NULL)
  .out
}

#' Build the 2nd-order augmented model (state + 1st/2nd-order sensitivities and
#' the prediction chain f1/f2).  This is the middle robustness tier: it omits the
#' O(neta^3) third-order states of `.foceiAnalyticAugModel`, so it is a much
#' smaller ODE that is far more likely to solve; the third-order tensor is then
#' recovered by finite-differencing the analytic 2nd-order sensitivities.
#' @return list(augMod, etav, neta, st, P2) or `NULL` on failure
#' @noRd
.foceiAnalyticSens2Model <- function(ui) {
  tryCatch({
    s <- ui$foceiEtaS; st <- rxode2::rxStateOde(s); neta <- s$..maxEta
    etav <- paste0("ETA_", seq_len(neta), "_")
    s1 <- s$..sens; s2 <- rxode2::.rxSens(s, etav, etav)
    .Dn <- function(e, v) symengine::D(e, symengine::S(v))
    .sn1 <- function(j, ...) symengine::S(paste0("rx__sens_", j, "_BY_", paste(c(...), collapse = "_BY_"), "__"))
    .toRx <- function(.l) rxode2::rxFromSE(.l)
    .g1 <- function(ex, p) { e <- .Dn(ex, p); for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p); e }
    .g2 <- function(ex, p, q) { gq <- .g1(ex, q); e <- .Dn(gq, p)
      for (k in st) e <- e + .Dn(gq, k) * .sn1(k, p); for (j in st) e <- e + .Dn(ex, j) * .sn1(j, p, q); e }
    pred <- get("rx_pred_", s); P2 <- expand.grid(i = etav, j = etav, stringsAsFactors = FALSE)
    fL1 <- vapply(etav, function(p) paste0("f1_", p, "=", .toRx(.g1(pred, p))), character(1))
    fL2 <- vapply(seq_len(nrow(P2)), function(r)
      paste0("f2_", P2$i[r], "_", P2$j[r], "=", .toRx(.g2(pred, P2$i[r], P2$j[r]))), character(1))
    baseOde <- vapply(st, function(x) paste0("d/dt(", x, ")=", .toRx(get(paste0("rx__d_dt_", x, "__"), s))), character(1))
    modTxt <- paste(c(baseOde, s1, s2, paste0("predf=", .toRx(pred)), fL1, fL2), collapse = "\n")
    modTxt <- gsub("ETA\\[([0-9]+)\\]", "ETA_\\1_", modTxt); modTxt <- gsub("THETA\\[([0-9]+)\\]", "THETA_\\1_", modTxt)
    list(augMod = rxode2::rxode2(modTxt), etav = etav, neta = neta, st = st, P2 = P2)
  }, error = function(e) NULL)
}

#' Harmonic-mean gate for the Shi (2021) adaptive step over a vector-valued
#' function (faithful to nlmixr2's `shi21.cpp`): one third-difference ratio per
#' component is collapsed to a single value by their harmonic mean, with a
#' correction for components that are locally linear (zero third difference).
#' @noRd
.shiHarmonicMean <- function(allv) {
  if (length(allv) == 1L) return(allv[1L])
  zero <- allv == 0; n <- sum(!zero); nzero <- sum(zero)
  if (n == 0L) return(0)
  correction <- (n - nzero) / n; if (correction <= 0) correction <- 1
  (n / sum(1 / allv[!zero])) * correction
}

#' One Shi (2021) central-difference probe at step `h` (paper Algorithm 3.1):
#' the harmonic-mean third-difference ratio plus the +/-h evaluations, or a
#' non-finite marker naming the evaluation that failed.
#' @noRd
.shiRC <- function(h, ffun, ef, t, idx, f0) {
  pert <- function(d) { tt <- t; tt[idx] <- tt[idx] + d; ffun(tt) }
  fp1 <- pert(h);      if (is.null(fp1)) return(list(r = NA_real_, bad = "p1"))
  fm1 <- pert(-h);     if (is.null(fm1)) return(list(r = NA_real_, bad = "m1", fp1 = fp1))
  fp3 <- pert(3 * h);  if (is.null(fp3)) return(list(r = NA_real_, bad = "p3", fp1 = fp1, fm1 = fm1))
  fm3 <- pert(-3 * h); if (is.null(fm3)) return(list(r = NA_real_, bad = "m3", fp1 = fp1, fm1 = fm1))
  allv <- abs(fp3 - 3 * fp1 + 3 * fm1 - fm3) / (8 * ef)
  list(r = .shiHarmonicMean(allv), fp1 = fp1, fm1 = fm1)
}

#' Shi (2021) adaptive central-difference step + gradient for a vector function
#' `ffun` along coordinate `idx` (R port of nlmixr2 `shi21Central`).  The step is
#' driven down/up by the third-difference ratio until it sits in [rl, ru];
#' returns the accepted step and the central-difference gradient at it.
#' @noRd
.shi21CentralR <- function(ffun, t, idx, f0, ef = 7e-7, rl = 1.5, ru = 4.5, nu = 3.0, maxiter = 15L) {
  h <- (3 * ef)^(1 / 3); l <- 0; u <- Inf; hlast <- h; gr <- NULL
  for (iter in seq_len(maxiter)) {
    sc <- .shiRC(h, ffun, ef, t, idx, f0)
    if (is.na(sc$r)) {                                   # an evaluation was non-finite
      if (sc$bad == "p1") { h <- h * 0.5 / 3; next }
      if (sc$bad == "m1") { if (is.null(gr)) gr <- (sc$fp1 - f0) / h; h <- h * 0.5 / 3; next }
      h <- h * 2 / 3; if (is.null(gr)) { gr <- (sc$fp1 - sc$fm1) / (2 * h); hlast <- h }; next
    }
    gr <- (sc$fp1 - sc$fm1) / (2 * h); hlast <- h
    if (sc$r < rl) l <- h else if (sc$r > ru) u <- h else break
    if (!is.finite(u)) h <- nu * h else if (l == 0) h <- h / nu else h <- (l + u) / 2
  }
  list(h = hlast, gr = gr)
}

#' Solve the 2nd-order model and obtain the 3rd-order tensor `Ath = d3f/deta3` by
#' SHI (2021) ADAPTIVE central differences of the analytic 2nd-order
#' sensitivities `A = d2f/deta2` (middle robustness tier).  The step is chosen per
#' eta-direction with the harmonic-mean gate over all components of `A`, exactly
#' as nlmixr2's inner-problem differences (`shi21Central` in `inner.cpp`) -- no
#' fixed step to misjudge on awkward/stiff models.  Then symmetrize (the 3rd
#' derivative is fully symmetric).  Returns the same `list(f,a,A,Ath)` the exact
#' path returns, so the assembly is identical; `NULL` on any failure.
#' @noRd
.foceiAnalyticSolveSubjectFD3 <- function(m2, params, ev, times, etav, neta, P2, ehat, ef = 7e-7) {
  base <- params[!names(params) %in% etav]
  solveA <- function(eta) {
    .d <- tryCatch(as.data.frame(rxode2::rxSolve(m2$augMod, params = c(base, setNames(eta, etav)), ev,
        returnType = "data.frame", atol = 1e-10, rtol = 1e-10)), error = function(e) NULL, warning = function(w) NULL)
    if (is.null(.d)) return(NULL); .d <- .d[.d$time %in% times, , drop = FALSE]
    if (nrow(.d) == 0L || !all(c("predf", paste0("f1_", etav)) %in% names(.d))) return(NULL)
    A <- array(0, c(nrow(.d), neta, neta))
    for (r in seq_len(nrow(P2))) A[, match(P2$i[r], etav), match(P2$j[r], etav)] <- .d[[paste0("f2_", P2$i[r], "_", P2$j[r])]]
    list(f = .d$predf, a = matrix(sapply(etav, function(p) .d[[paste0("f1_", p)]]), nrow(.d), neta), A = A)
  }
  Aflat <- function(eta) { E <- solveA(eta); if (is.null(E) || !all(is.finite(E$A))) return(NULL); as.vector(E$A) }
  E0 <- solveA(ehat); if (is.null(E0)) return(NULL); nobs <- length(E0$f); f0 <- as.vector(E0$A)
  Ath <- array(0, c(nobs, neta, neta, neta))
  for (n in seq_len(neta)) {
    sc <- .shi21CentralR(Aflat, ehat, n, f0, ef = ef)
    if (is.null(sc$gr) || !all(is.finite(sc$gr))) return(NULL)
    Ath[, , , n] <- array(sc$gr, c(nobs, neta, neta))
  }
  S <- array(0, c(nobs, neta, neta, neta))
  pr <- list(c(1,2,3), c(1,3,2), c(2,1,3), c(2,3,1), c(3,1,2), c(3,2,1))
  for (l in seq_len(neta)) for (m in seq_len(neta)) for (n in seq_len(neta)) {
    ix <- c(l, m, n); S[, l, m, n] <- Reduce(`+`, lapply(pr, function(p) Ath[, ix[p[1]], ix[p[2]], ix[p[3]]])) / 6
  }
  .out <- list(f = E0$f, a = E0$a, A = E0$A, Ath = S)
  if (!all(is.finite(.out$f)) || !all(is.finite(.out$a)) || !all(is.finite(.out$A)) || !all(is.finite(.out$Ath))) return(NULL)
  .out
}

#' Per-subject FULL observed-information R-matrix block over all population
#' parameters: structural (mu-ref) theta, residual sigma, and Omega (diagonal or
#' block) via the `omd` non-Cholesky variance-covariance derivatives.
#'
#' Pure assembly from already-evaluated sensitivities (`E$a/A/Ath`, the residual
#' params via `ef`); ported from the validated prototype `inst/tools/rmatrix-full.R`.
#' Each entry is the data term (envelope/Schur) plus the log-determinant term
#' (moving mode, 3rd-order sensitivities), generalized to the parameter types via
#' the `typ()`/`Mcol()`/`Smat()`/... family.  Parameter order: `nth` structural
#' thetas, `nsg` sigma (`sgVar`), then `neta` Omega variances.  Diagonal Omega.
#' @noRd
.foceiAnalyticSubjectR <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, omd, covFull = FALSE,
                                   ndir = neta, dirTh = seq_len(nth)) {
  # Observed-information R for one subject.  The inner Hessian H is over the etas
  # ("eta-slots", 1:neta); each structural theta differentiates in its own
  # "direction-slot" (1:ndir): a mu-ref theta reuses its eta direction, an eta-less
  # theta gets its own added direction (dirTh maps theta index -> direction).  For a
  # fully mu-ref model ndir == neta, dirTh == 1:nth, and this reduces to the original.
  tr <- function(M) sum(diag(M))
  a <- E$a; A <- E$A; Ath <- E$Ath; f <- E$f; y <- E$y; Oi <- solve(Om)
  evf <- function(e) ef$ev(e, f, y)
  rd <- list(r1 = evf(ef$sc$r1), r2 = evf(ef$sc$r2), r3 = evf(ef$sc$r3))
  pf <- list(p = evf(ef$sc$p), p1 = evf(ef$sc$p1), p2 = evf(ef$sc$p2))
  np <- nth + nsg + omd$nom
  ei <- seq_len(neta); di <- seq_len(ndir); ae <- a[, ei, drop = FALSE]  # eta-cols for sigma/Omega
  dir <- function(p) dirTh[p]
  H <- Oi; for (l in ei) for (m in ei) H[l,m] <- H[l,m] + sum(rd$r2*a[,l]*a[,m] + rd$r1*A[,l,m])
  N <- matrix(0, neta, ndir); for (l in ei) for (d in di) N[l,d] <- sum(rd$r2*a[,l]*a[,d] + rd$r1*A[,l,d])
  HiM <- solve(H)
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l,m] <- Ht[l,m] + sum(pf$p*a[,l]*a[,m]); Hti <- solve(Ht)
  ouAA <- function(v) { M <- matrix(0,neta,neta); for (l in ei) for (m in ei) M[l,m] <- sum(v*a[,l]*a[,m]); M }
  dHtD <- lapply(di, function(s) { D <- matrix(0,neta,neta); for (l in ei) for (m in ei)
    D[l,m] <- sum(pf$p1*a[,s]*a[,l]*a[,m] + pf$p*A[,l,s]*a[,m] + pf$p*a[,l]*A[,m,s]); D })
  d2HtDD <- lapply(di, function(s) lapply(di, function(t) { D <- matrix(0,neta,neta); for (l in ei) for (m in ei)
    D[l,m] <- sum(pf$p2*a[,s]*a[,t]*a[,l]*a[,m] +
      pf$p1*(A[,s,t]*a[,l]*a[,m]+a[,s]*A[,l,t]*a[,m]+a[,s]*a[,l]*A[,m,t]+a[,t]*A[,l,s]*a[,m]+a[,t]*a[,l]*A[,m,s]) +
      pf$p*(Ath[,l,s,t]*a[,m]+A[,l,s]*A[,m,t]+A[,l,t]*A[,m,s]+a[,l]*Ath[,m,s,t])); D }))
  Cen <- vapply(ei, function(l) 0.5*tr(Hti %*% dHtD[[l]]), numeric(1))
  Cee <- matrix(0,neta,neta); for (s in ei) for (t in ei)
    Cee[s,t] <- 0.5*(tr(Hti %*% d2HtDD[[s]][[t]]) - tr(Hti %*% dHtD[[s]] %*% Hti %*% dHtD[[t]]))
  Tn <- array(0,c(neta,ndir,ndir)); for (l in ei) for (s in di) for (t in di)
    Tn[l,s,t] <- sum(rd$r3*a[,l]*a[,s]*a[,t] + rd$r2*(A[,l,s]*a[,t]+A[,l,t]*a[,s]+A[,s,t]*a[,l]) + rd$r1*Ath[,l,s,t])
  Ndd <- function(da,db) sum(rd$r2*a[,da]*a[,db] + rd$r1*A[,da,db])
  typ <- function(p) if (p <= nth) "th" else if (p <= nth+nsg) "sg" else "om"
  sgi <- function(p) sgVar[p - nth]
  omc <- function(p) p - nth - nsg
  PVper <- function(p) lapply(ef$per[[sgi(p)]], evf)
  PVpair <- function(aa, bb) { s1 <- sgi(aa); s2 <- sgi(bb); key <- if (paste0(s1,s2) %in% names(ef$pair)) paste0(s1,s2) else paste0(s2,s1)
    lapply(ef$pair[[key]], evf) }
  # Omega enters only the prior (1/2 eta' Omega^-1 eta + 1/2 ln|Omega|) and H~'s
  # +Omega^-1 term, so every Omega-block quantity is an E-basis contraction from
  # `omd` (non-Cholesky variance-covariance derivatives) -- diagonal or block.
  Mcol <- function(p) { t <- typ(p)
    if (t=="th") return(N[,dir(p)]); if (t=="sg") return(as.numeric(crossprod(ae, PVper(p)$rf)))
    as.numeric(omd$dOi[[omc(p)]] %*% ehat) }
  dHt_p <- function(p) { t <- typ(p)
    if (t=="th") return(dHtD[[dir(p)]]); if (t=="sg") return(ouAA(PVper(p)$ps))
    omd$dOi[[omc(p)]] }
  d2HtEtaP <- function(p,l) { t <- typ(p)
    if (t=="th") return(d2HtDD[[dir(p)]][[l]]); if (t=="om") return(matrix(0,neta,neta))
    P <- PVper(p); D <- matrix(0,neta,neta); for (s in ei) for (m in ei)
      D[s,m] <- sum(P$pf*a[,l]*a[,s]*a[,m] + P$ps*(A[,s,l]*a[,m]+a[,s]*A[,m,l])); D }
  d2Ht_pp <- function(aa,bb) { ta<-typ(aa); tb<-typ(bb)
    if (ta=="th"&&tb=="th") return(d2HtDD[[dir(aa)]][[dir(bb)]])
    if (ta=="om"&&tb=="om") return(omd$d2Oi[[omc(aa)]][[omc(bb)]])
    if (ta=="om"||tb=="om") return(matrix(0,neta,neta))
    if (ta=="sg"&&tb=="sg") return(ouAA(PVpair(aa,bb)$pss))
    thp<-if(ta=="th")aa else bb; sg<-if(ta=="th")bb else aa; d2HtEtaP(sg, dir(thp)) }
  Smat <- function(p) { t<-typ(p)
    if (t=="th") { M<-matrix(0,neta,ndir); for (l in ei) for (s in di) M[l,s]<-Tn[l,dir(p),s]; return(M) }
    if (t=="om") return(omd$dOi[[omc(p)]])
    P<-PVper(p); M<-matrix(0,neta,ndir); for (l in ei) for (s in di) M[l,s]<-sum(P$rff*a[,s]*a[,l]+P$rf*A[,l,s]); M }
  Svec <- function(aa,bb) { ta<-typ(aa); tb<-typ(bb)
    if (ta=="th"&&tb=="th") return(Tn[,dir(aa),dir(bb)])
    if (ta=="om"&&tb=="om") return(as.numeric(omd$d2Oi[[omc(aa)]][[omc(bb)]] %*% ehat))
    if (ta=="om"||tb=="om") return(rep(0,neta))
    if (ta=="sg"&&tb=="sg") return(as.numeric(crossprod(ae, PVpair(aa,bb)$rfss)))
    thp<-if(ta=="th")aa else bb; sg<-if(ta=="th")bb else aa; Smat(sg)[,dir(thp)] }
  d2Phi <- function(aa,bb) { ta<-typ(aa); tb<-typ(bb)
    if (ta=="th"&&tb=="th") return(Ndd(dir(aa),dir(bb)))
    if (ta=="om"&&tb=="om") return(0.5*as.numeric(t(ehat)%*%omd$d2Oi[[omc(aa)]][[omc(bb)]]%*%ehat) + 0.5*omd$d2LD[omc(aa),omc(bb)])
    if (ta=="om"||tb=="om") return(0)
    if (ta=="sg"&&tb=="sg") return(sum(PVpair(aa,bb)$rss))
    thp<-if(ta=="th")aa else bb; sg<-if(ta=="th")bb else aa; as.numeric(crossprod(a[,dir(thp)], PVper(sg)$rf)) }
  # covFull=FALSE (theta-only, the default scope): the theta-theta block still needs
  # Tn / d2HtDD (3rd-order sensitivities), so the augmented ODE solve is unchanged;
  # only the sigma/Omega columns of the assembly are skipped (np -> nth).
  npE <- if (covFull) np else nth
  etaP <- matrix(vapply(1:npE, function(p) as.numeric(-HiM %*% Mcol(p)), numeric(neta)), nrow = neta)  # neta x npE (neta==1 safe)
  eta2 <- function(aa,bb) { b <- Svec(aa,bb) + Smat(aa)[,ei,drop=FALSE]%*%etaP[,bb] + Smat(bb)[,ei,drop=FALSE]%*%etaP[,aa]
    for (l in ei) b[l] <- b[l] + as.numeric(t(etaP[,aa])%*%Tn[l,ei,ei]%*%etaP[,bb]); as.numeric(-HiM%*%b) }
  Cpe <- function(p,l) 0.5*(tr(Hti%*%d2HtEtaP(p,l)) - tr(Hti%*%dHt_p(p)%*%Hti%*%dHtD[[l]]))
  Cpp <- function(aa,bb) 0.5*(tr(Hti%*%d2Ht_pp(aa,bb)) - tr(Hti%*%dHt_p(aa)%*%Hti%*%dHt_p(bb)))
  R <- matrix(0,npE,npE)
  for (aa in 1:npE) for (bb in aa:npE) {              # R is symmetric -- fill upper, mirror
    dat <- d2Phi(aa,bb) - as.numeric(t(Mcol(aa))%*%HiM%*%Mcol(bb))
    ld <- Cpp(aa,bb) + sum(vapply(ei,function(l)Cpe(aa,l),numeric(1))*etaP[,bb]) +
          sum(vapply(ei,function(mm)Cpe(bb,mm),numeric(1))*etaP[,aa]) +
          as.numeric(t(etaP[,aa])%*%Cee%*%etaP[,bb]) + sum(Cen*eta2(aa,bb))
    R[aa,bb] <- R[bb,aa] <- dat + ld
  }
  R
}

#' Full analytic FOCEI/FOCE covariance for a fitted nlmixr2 object.
#'
#' Computes the exact observed-information R-matrix (data term + log-determinant
#' term) over ALL population parameters -- the structural (mu-referenced) fixed
#' effects, the residual sigma, and the `Omega` variances AND covariances
#' (diagonal or block, via the non-Cholesky `rxOmegaVarCovDeriv` derivatives) --
#' and returns the covariance `solve(R)`.  Finite-difference-free; runs on
#' released rxode2 (pure-R fallbacks are used when `rxExpandSens3_` /
#' `rxOmegaVarCovDeriv` are absent).  Returns `NULL` (caller should fall back) if
#' the model is out of scope (non-mu-ref eta, unsupported error model) or any
#' per-subject augmented solve fails.
#'
#' @param fit a fitted nlmixr2 focei object
#' @param sens sensitivity source for the log-determinant term: `"exact3"`
#'   (analytic 3rd-order, the default) or `"fd2"` (analytic 2nd-order with finite
#'   differences of the 2nd-order sensitivities for the 3rd-order tensor -- a
#'   lighter ODE that solves more reliably, the middle robustness tier).  Both
#'   feed the identical R-matrix assembly.
#' @return list(cov, se, R, params), with `params` covering theta, sigma, the
#'   `om.<theta>` Omega variances and any `cov.<theta>.<theta>` covariances, or
#'   `NULL`
#' @noRd
foceiCovAnalytic <- function(fit, sens = c("exact3", "fd2"), covFull = FALSE) {
  sens <- match.arg(sens)
  ui <- fit$finalUi
  if (!.hasRxExpandSens3()) return(NULL)            # need rxExpandSens2_ + symengine
  ef <- .foceiAnalyticErrFull(ui)
  if (is.null(ef)) return(NULL)                     # unsupported error model -> fallback

  ini <- ui$iniDf
  muRef <- ui$muRefDataFrame                       # theta <-> eta pairing
  etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
  etaRows <- etaRows[order(etaRows$neta1), , drop = FALSE]
  etaNames <- etaRows$name
  neta <- length(etaNames)
  if (neta == 0L) return(NULL)
  thetaForEta <- muRef$theta[match(etaNames, muRef$eta)]
  if (anyNA(thetaForEta)) return(NULL)             # non-mu-ref eta -> out of scope
  fx <- function(nm) { if (is.null(ini$fix)) return(rep(FALSE, length(nm)))
    i <- match(nm, ini$name); !is.na(i) & !is.na(ini$fix[i]) & ini$fix[i] }
  if (any(fx(thetaForEta))) return(NULL)           # fixed structural theta breaks eta indexing -> builtin
  keep <- !fx(ef$sgName); ef$sgVar <- ef$sgVar[keep]; ef$sgName <- ef$sgName[keep]  # drop fixed sigma
  Om <- fit$omega
  blocks <- .omegaBlocks(Om)                       # free Omega lower-triangle (diagonal or block)
  pairs <- do.call(rbind, lapply(blocks, function(b)
    do.call(rbind, lapply(seq_along(b), function(a) cbind(b[a], b[seq_len(a)])))))
  pairs <- pairs[!.omegaFixed(ini, pairs), , drop = FALSE]   # drop fixed Omega elements
  omd <- .omegaVarCovDeriv(Om, pairs)
  thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]
  thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
  # Uniform direction assembly: every ESTIMATED structural theta is handled as ONE
  # direction.  A mu-ref theta reuses its eta's direction (free -- d/dtheta == d/deta).
  # Any other structural theta -- a covariate coefficient OR an eta-less theta (a `tv`
  # with no eta, a transit `mtt`, ...) -- gets its OWN THETA_j_ direction carrying its
  # TRUE sensitivity.  Nothing is "detected" or scaled: correctness is by construction,
  # so there is no covariate heuristic to be fragile.  (sigma + Omega as before.)
  thStructRows <- thRows[!fx(thRows$name) & !(thRows$name %in% ef$sgName), , drop = FALSE]
  thStruct <- thStructRows$name
  nth <- length(thStruct)
  if (nth == 0L) return(NULL)
  etaDirs <- paste0("ETA_", seq_len(neta), "_")
  dirTh <- integer(nth); nonMuTheta <- character(0)
  for (p in seq_len(nth)) {
    k <- match(thStruct[p], thetaForEta)             # eta this theta mu-references (NA if none)
    if (!is.na(k)) {
      dirTh[p] <- k                                  # reuse eta k's direction (free)
    } else {
      nonMuTheta <- c(nonMuTheta, paste0("THETA_", thStructRows$ntheta[p], "_"))
      dirTh[p] <- neta + length(nonMuTheta)          # its own true-sensitivity direction
    }
  }
  dirs <- c(etaDirs, nonMuTheta)
  ndir <- length(dirs)
  # fd2 tier is only wired for the pure mu-ref case (dirs = etas); a non-mu-ref theta
  # needs the exact3 direction-set model -> bow out of fd2 so the ladder uses exact3
  if (sens == "fd2" && length(nonMuTheta) > 0L) return(NULL)

  th <- setNames(thRows$est, paste0("THETA_", seq_len(nrow(thRows)), "_"))
  ebes <- as.matrix(fit$eta[, etaNames, drop = FALSE])
  nsg <- length(ef$sgVar); np <- nth + nsg + omd$nom

  am <- if (sens == "fd2") .foceiAnalyticSens2Model(ui) else .foceiAnalyticAugModelDirs(ui, dirs)
  if (is.null(am)) return(NULL)
  if (sens == "fd2") { if (am$neta != neta) return(NULL) } else if (am$ndir != ndir) return(NULL)

  ds <- fit$dataSav
  ids <- fit$eta$ID
  npOut <- if (covFull) np else nth                  # covFull=FALSE: structural thetas only

  R <- matrix(0, npOut, npOut)
  for (i in seq_along(ids)) {
    s <- ds[ds$ID == ids[i], , drop = FALSE]         # solve over the subject's ACTUAL events
    obs <- s[s$EVID == 0, , drop = FALSE]             # (CMT/EVID/II/SS/ADDL/covariates preserved)
    p <- c(th, setNames(ebes[i, ], etaDirs))
    E <- if (sens == "fd2")
      .foceiAnalyticSolveSubjectFD3(am, p, s, obs$TIME, etaDirs, neta, am$P2, ebes[i, ])
    else
      .foceiAnalyticSolveDir(am, p, s, obs$TIME)
    if (is.null(E)) return(NULL)                   # solve failure -> next tier (caller)
    E$y <- obs$DV
    Ri <- tryCatch(.foceiAnalyticSubjectR(E, ebes[i, ], Om, ef, neta, nth, nsg, ef$sgVar, omd,
                                          covFull = covFull, ndir = ndir, dirTh = dirTh),
                   error = function(e) NULL)
    if (is.null(Ri) || !all(is.finite(Ri))) return(NULL)
    R <- R + Ri
  }
  cov <- tryCatch(solve(R), error = function(e) NULL)
  if (is.null(cov)) return(NULL)
  omNm <- apply(pairs, 1, function(p) if (p[1] == p[2]) paste0("om.", thetaForEta[p[1]])
                else paste0("cov.", thetaForEta[p[1]], ".", thetaForEta[p[2]]))
  nm <- if (covFull) c(thStruct, ef$sgName, omNm) else thStruct
  dimnames(R) <- dimnames(cov) <- list(nm, nm)
  list(cov = cov, se = setNames(suppressWarnings(sqrt(diag(cov))), nm), R = R, params = nm)  # NaN flags non-PD
}

