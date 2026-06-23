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
# falls back to the existing finite-difference covariance.  The analytic path is
# therefore never partially applied.

#' Pure-R generator of the 3rd-order forward sensitivity symengine code.
#'
#' Faithful port of rxode2's C++ `rxExpandSens3_` (string generation only, no
#' ODE/SUNDIALS dependency), so the analytic covariance runs on a released
#' rxode2 (>= 5.0.2 has `rxExpandSens2_` + symengine) without the third-order
#' C++ function.  When rxode2 ships `rxExpandSens3_` (PR nlmixr2/rxode2#1092)
#' that faster version is used instead (see [.rxExpandSens3]).
#'
#' Builds, for each (state cS, vars pa,pb,pc), the total d/d(pa) of the stored
#' 2nd-order sensitivity RHS `v2 = d/dt(s^{pb,pc})`, with the unique{pb,pc}
#' dedup that avoids a spurious x2 when pb == pc.
#' @return data.frame with columns ddt, ddtS, ddS2, line (matches rxExpandSens3_)
#' @noRd
.rxExpandSens3R <- function(state, s1, s2, s3) {
  .res <- function(v) if (v %in% c("e","E","EulerGamma","Catalan","GoldenRatio","I")) paste0("rx_SymPy_Res_", v) else v
  g <- expand.grid(cS = state, pa = s1, pb = s2, pc = s3, stringsAsFactors = FALSE)
  nm <- function(cS, ...) paste0("rx__sens_", cS, "_BY_", paste(c(...), collapse = "_BY_"), "__")
  ddtn <- function(cS, ...) paste0("rx__d_dt_", nm(cS, ...), "__")
  line <- character(nrow(g)); ddt <- character(nrow(g)); ddtS <- character(nrow(g)); ddS2 <- character(nrow(g))
  for (r in seq_len(nrow(g))) {
    cS <- g$cS[r]; pa <- g$pa[r]; pb <- g$pb[r]; pc <- g$pc[r]
    sp <- nm(cS, pa, pb, pc); ddt[r] <- paste0("d/dt(", sp, ")"); ddS2[r] <- sp; ddtS[r] <- ddtn(cS, pa, pb, pc)
    v2 <- ddtn(cS, pb, pc)
    L <- paste0("D(", v2, ",\"", .res(pa), "\")")
    for (xm in state)
      L <- paste0(L, "+D(", v2, ",\"", .res(xm), "\")*", nm(xm, pa),
                  "+rx__df_", cS, "_dy_", xm, "__*", nm(xm, pa, pb, pc))
    for (xm in state)
      L <- paste0(L, "+D(", v2, ",\"", nm(xm, pb), "\")*", nm(xm, pa, pb))
    if (pc != pb) for (xm in state)
      L <- paste0(L, "+D(", v2, ",\"", nm(xm, pc), "\")*", nm(xm, pa, pc))
    line[r] <- paste0("assign(\"", ddtS[r], "\", with(model,", L, "), envir=model)")
  }
  data.frame(ddt = ddt, ddtS = ddtS, ddS2 = ddS2, line = line, stringsAsFactors = FALSE)
}

#' Third-order sensitivity generator: rxode2's C++ version if present, else the
#' pure-R port. The analytic covariance therefore needs no special rxode2 build.
#' @noRd
.rxExpandSens3 <- function(state, s1, s2, s3) {
  if (exists("rxExpandSens3_", envir = asNamespace("rxode2"), inherits = FALSE)) {
    return(rxode2::rxExpandSens3_(state, s1, s2, s3))
  }
  .rxExpandSens3R(state, s1, s2, s3)
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

#' Decide whether a model is in-scope for the analytic covariance Hessian.
#'
#' In scope (initial release): Gaussian add/prop/combined residual error,
#' mu-referenced structural thetas, and a model rxode2 can build 3rd-order
#' sensitivities for.  Out of scope (-> FD fallback): generalized likelihood,
#' BLOQ/censoring (M3/M4), residual-error random effects, IOV, etc.
#' @param ui rxode2 UI object
#' @return `TRUE` if in scope, otherwise `FALSE`
#' @noRd
.foceiAnalyticInScope <- function(ui) {
  if (!.hasRxExpandSens3()) return(FALSE)
  tryCatch({
    # in scope only for a supported Gaussian add/prop/combined error model
    if (is.null(.foceiAnalyticErrModel(ui))) return(FALSE)
    # mu-referenced structural thetas required (every eta paired with a theta)
    ini <- ui$iniDf
    etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
    if (nrow(etaRows) == 0L) return(FALSE)
    if (anyNA(ui$muRefDataFrame$theta[match(etaRows$name, ui$muRefDataFrame$eta)])) return(FALSE)
    TRUE   # both FOCE and FOCEI interaction modes are in scope
  }, error = function(e) FALSE)
}

#' Build the error-model scalar closures rho(f,y) and p(f) and their exact
#' f-derivatives, from the residual variance R(f) (a quoted expression in `f`
#' plus named constants).  rho = 1/2((y-f)^2/R + log R); p = 1/R + 1/2 (R'/R)^2
#' (the FOCEI approximate-Hessian weight).  Generic over add / prop / combined.
#' @noRd
.foceiAnalyticErrScalars <- function(Rexpr, consts) {
  rhoE <- substitute(0.5 * ((y - f)^2 / RR + log(RR)), list(RR = Rexpr))
  dR   <- D(Rexpr, "f")
  pE   <- substitute(1 / RR + 0.5 * (DR / RR)^2, list(RR = Rexpr, DR = dR))
  r1 <- D(rhoE, "f"); r2 <- D(r1, "f"); r3 <- D(r2, "f")
  p1 <- D(pE, "f");  p2 <- D(p1, "f")
  ev <- function(e, f, y) eval(e, c(list(f = f, y = y), as.list(consts)))
  list(
    rho = function(f, y) list(r1 = ev(r1, f, y), r2 = ev(r2, f, y), r3 = ev(r3, f, y)),
    p   = function(f)    list(p = ev(pE, f, 0), p1 = ev(p1, f, 0), p2 = ev(p2, f, 0))
  )
}

#' Residual-variance expression + constants from a fitted nlmixr2 UI.
#' Supports additive, proportional, and combined (combined2 / combined1) Gaussian
#' error.  Returns `NULL` for out-of-scope error models (-> FD fallback).
#' @noRd
.foceiAnalyticErrModel <- function(ui) {
  ini <- ui$iniDf
  errRows <- ini[!is.na(ini$err), , drop = FALSE]
  add <- errRows$est[errRows$err %in% c("add", "lnorm")]
  prop <- errRows$est[errRows$err %in% c("prop", "propT", "propF")]
  addPr <- tryCatch(rxode2::rxGetControl(ui, "addProp", "combined2"), error = function(e) "combined2")
  if (length(add) == 1L && length(prop) == 0L)
    return(list(R = quote(a2), consts = c(a2 = add^2)))
  if (length(add) == 0L && length(prop) == 1L)
    return(list(R = quote(p2 * f^2), consts = c(p2 = prop^2)))
  if (length(add) == 1L && length(prop) == 1L) {
    if (identical(addPr, "combined1"))                       # sd = add + prop*f
      return(list(R = quote((sqrt(a2) + sqrt(p2) * f)^2), consts = c(a2 = add^2, p2 = prop^2)))
    return(list(R = quote(a2 + p2 * f^2), consts = c(a2 = add^2, p2 = prop^2)))  # combined2
  }
  NULL
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
    model <- s                                  # rxExpandSens3_ lines reference `model`
    st   <- rxode2::rxStateOde(s)
    neta <- s$..maxEta
    etav <- paste0("ETA_", seq_len(neta), "_")
    s1 <- s$..sens
    s2 <- rxode2::.rxSens(s, etav, etav)
    grd3 <- .rxExpandSens3(st, etav, etav, etav)
    invisible(lapply(c(grd3$ddtS, grd3$ddS2),
                     function(x) assign(x, symengine::Symbol(x), envir = s)))
    s3 <- vapply(seq_len(nrow(grd3)), function(i) {
      .l <- eval(parse(text = grd3$line[i]))
      paste0(grd3$ddt[i], "=", rxode2::rxFromSE(.l))
    }, character(1))

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

#' Per-subject observed-information R-matrix block (data term + log-determinant
#' term).
#'
#' Pure assembly from already-evaluated sensitivities; ported verbatim from the
#' validated prototype (`inst/tools/rmatrix-full.R`).  Writing the FOCE/FOCEI
#' objective as a data term `l_i(etaHat)` plus the log-determinant of the
#' first-order inner Hessian, the population Hessian splits into the data-term
#' contribution (reusing the 2nd-order sensitivities of the gradient) and the
#' log-determinant-term contribution (which carries the 3rd-order sensitivities).
#' `Om` is the random-effects covariance, `ehat` the subject EBE, `th` the
#' residual params.
#' @noRd
.foceiAnalyticSubjectR <- function(ev, ehat, th, Om, sigPart, neta) {
  .tr <- function(M) sum(diag(M))
  a <- ev$a; A <- ev$A; Ath <- ev$Ath; f <- ev$f
  Oi <- solve(Om)
  rd <- sigPart$rho(f)              # list r1,r2,r3 : rho',rho'',rho''' per obs
  pf <- sigPart$p(f)               # list p,p1,p2 : FOCEI approx-Hessian weight
  # true inner Hessian H and N = H - Omega^{-1}
  H <- Oi
  for (l in 1:neta) for (m in 1:neta)
    H[l, m] <- H[l, m] + sum(rd$r2 * a[, l] * a[, m] + rd$r1 * A[, l, m])
  N <- H - Oi; HiM <- solve(H)
  # data term: profile Hessian of l_i (envelope/Schur complement)
  R_data <- N - N %*% solve(H) %*% N
  # H~ and its eta-derivatives (log-determinant-term machinery)
  Ht <- Oi
  for (l in 1:neta) for (m in 1:neta) Ht[l, m] <- Ht[l, m] + sum(pf$p * a[, l] * a[, m])
  dHt <- lapply(1:neta, function(s) { D <- matrix(0, neta, neta)
    for (l in 1:neta) for (m in 1:neta)
      D[l, m] <- sum(pf$p1 * a[, s] * a[, l] * a[, m] + pf$p * A[, l, s] * a[, m] + pf$p * a[, l] * A[, m, s]); D })
  d2Ht <- lapply(1:neta, function(s) lapply(1:neta, function(t) { D <- matrix(0, neta, neta)
    for (l in 1:neta) for (m in 1:neta)
      D[l, m] <- sum(pf$p2 * a[, s] * a[, t] * a[, l] * a[, m] +
        pf$p1 * (A[, s, t] * a[, l] * a[, m] + a[, s] * A[, l, t] * a[, m] + a[, s] * a[, l] * A[, m, t] +
                 a[, t] * A[, l, s] * a[, m] + a[, t] * a[, l] * A[, m, s]) +
        pf$p * (Ath[, l, s, t] * a[, m] + A[, l, s] * A[, m, t] + A[, l, t] * A[, m, s] + a[, l] * Ath[, m, s, t])); D }))
  Hti <- solve(Ht)
  Cen <- vapply(1:neta, function(l) 0.5 * .tr(Hti %*% dHt[[l]]), numeric(1))
  Cee <- matrix(0, neta, neta)
  for (s in 1:neta) for (t in 1:neta)
    Cee[s, t] <- 0.5 * (.tr(Hti %*% d2Ht[[s]][[t]]) - .tr(Hti %*% dHt[[s]] %*% Hti %*% dHt[[t]]))
  Tn <- array(0, c(neta, neta, neta))
  for (l in 1:neta) for (m in 1:neta) for (n in 1:neta)
    Tn[l, m, n] <- sum(rd$r3 * a[, l] * a[, m] * a[, n] +
      rd$r2 * (A[, l, m] * a[, n] + A[, l, n] * a[, m] + A[, m, n] * a[, l]) + rd$r1 * Ath[, l, m, n])
  etaTh <- -solve(H, N)
  mode2 <- array(0, c(neta, neta, neta))
  for (j in 1:neta) for (k in 1:neta) {
    b <- Tn[, j, k] + Tn[, , j] %*% etaTh[, k] + Tn[, , k] %*% etaTh[, j]
    for (l in 1:neta) b[l] <- b[l] + as.numeric(t(etaTh[, j]) %*% Tn[l, , ] %*% etaTh[, k])
    mode2[, j, k] <- -solve(H, b)
  }
  # log-determinant term: carried through the moving mode (3rd-order sensitivities)
  R_logDet <- matrix(0, neta, neta)
  for (j in 1:neta) for (k in 1:neta)
    R_logDet[j, k] <- Cee[j, k] + sum(Cee[j, ] * etaTh[, k]) + sum(Cee[k, ] * etaTh[, j]) +
      as.numeric(t(etaTh[, j]) %*% Cee %*% etaTh[, k]) + sum(Cen * mode2[, j, k])
  # NB: structural-theta block shown; sigma/Omega blocks add via the general
  # parameter assembly (inst/tools/rmatrix-full.R) -- wired in the full builder.
  R_data + R_logDet
}

#' Analytic FOCEI covariance (structural-theta block) for a fitted nlmixr2 object.
#'
#' Computes the exact observed-information R-matrix (data term + log-determinant
#' term) for the
#' structural (mu-referenced) fixed effects, finite-difference-free, and returns
#' the covariance `solve(R)`.  Runs on released rxode2 (uses the pure-R 3rd-order
#' generator when `rxExpandSens3_` is absent).  Returns `NULL` (caller should
#' fall back to the finite-difference covariance) if the model is out of scope or
#' any per-subject augmented solve fails.
#'
#' @param fit a fitted nlmixr2 focei object
#' @return list(cov, se, R, params) labelled by structural theta name, or `NULL`
#' @export
foceiCovAnalytic <- function(fit) {
  ui <- fit$finalUi
  if (!.foceiAnalyticInScope(ui)) return(NULL)
  em <- .foceiAnalyticErrModel(ui)
  if (is.null(em)) return(NULL)
  sig <- .foceiAnalyticErrScalars(em$R, em$consts)

  ini <- ui$iniDf
  muRef <- ui$muRefDataFrame                       # theta <-> eta pairing
  etaRows <- ini[!is.na(ini$neta1) & ini$neta1 == ini$neta2, , drop = FALSE]
  etaRows <- etaRows[order(etaRows$neta1), , drop = FALSE]
  etaNames <- etaRows$name
  neta <- length(etaNames)
  if (neta == 0L) return(NULL)
  # structural theta paired with each eta, in eta order
  thetaForEta <- muRef$theta[match(etaNames, muRef$eta)]
  if (anyNA(thetaForEta)) return(NULL)             # non-mu-ref eta -> out of scope
  # full theta value vector THETA_1.._n in ntheta order (incl error thetas)
  thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]
  thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
  th <- setNames(thRows$est, paste0("THETA_", seq_len(nrow(thRows)), "_"))
  Om <- fit$omega
  etav <- paste0("ETA_", seq_len(neta), "_")
  ebes <- as.matrix(fit$eta[, etaNames, drop = FALSE])

  am <- .foceiAnalyticAugModel(ui)
  if (is.null(am) || am$neta != neta) return(NULL)

  ds <- fit$dataSav
  ids <- fit$eta$ID
  R <- matrix(0, neta, neta)
  for (i in seq_along(ids)) {
    s <- ds[ds$ID == ids[i], , drop = FALSE]
    dose <- s[s$EVID != 0, , drop = FALSE]; obs <- s[s$EVID == 0, , drop = FALSE]
    ev <- rxode2::et()
    for (k in seq_len(nrow(dose))) ev <- ev |> rxode2::et(amt = dose$AMT[k], cmt = "depot", time = dose$TIME[k])
    ev <- ev |> rxode2::et(obs$TIME)
    E <- .foceiAnalyticSolveSubject(am$augMod, c(th, setNames(ebes[i, ], etav)),
                                    ev, obs$TIME, etav, neta, am$P2, am$P3)
    if (is.null(E)) return(NULL)                   # solve failure -> whole-population fallback
    yi <- obs$DV
    sp2 <- list(rho = function(f) sig$rho(f, yi), p = sig$p)
    Ri <- tryCatch(.foceiAnalyticSubjectR(E, ebes[i, ], th, Om, sp2, neta),
                   error = function(e) NULL)
    if (is.null(Ri) || !all(is.finite(Ri))) return(NULL)
    R <- R + Ri
  }
  cov <- tryCatch(solve(R), error = function(e) NULL)
  if (is.null(cov)) return(NULL)
  dimnames(R) <- dimnames(cov) <- list(thetaForEta, thetaForEta)
  list(cov = cov, se = setNames(sqrt(diag(cov)), thetaForEta), R = R, params = thetaForEta)
}

