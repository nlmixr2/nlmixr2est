# Analytic FOCE/FOCEI outer (population) gradient from Almquist (2015) sensitivity
# equations -- the first-derivative precursor of the analytic observed-information
# R-matrix in foceiCovAnalytic.R.  Enabled by foceiControl(fast=TRUE); called from
# the C++ outer-gradient hook (numericGrad -> analyticOuterGrad).  Returns a
# length-npars gradient of the objective (-2*logLik) in the model-theta scale, or
# NULL to fall back to the finite-difference gradient.
#
# The outer gradient needs at most SECOND-order state sensitivities (Almquist Eqs
# 38-40) -- one order less than the covariance R-matrix -- so it reuses the
# direction set, error machinery, and per-subject H/Ht/N/dHtD/etaP assembly from
# foceiCovAnalytic.R without the 3rd-order (Ath) tier.

#' Per-subject first-derivative (outer-gradient) contribution for FOCEI
#' (interaction=1).  Objective per subject (up to constants) is
#' `OFV_i = 2*Phi_i + log|Ht_i|`, where `Phi_i = -l_i*` is the data+prior part at
#' the EBE eta* and `Ht_i` is the Laplace (Eq-14) determinant Hessian.  By the
#' envelope theorem (`dl/deta = 0` at eta*) the data term is the explicit partial;
#' the log-determinant term carries the moving-mode `etaP = d eta*/d p` (Eq 46).
#'
#' Reuses the cov machinery's H/N/HiM/Ht/dHtD/etaP/Mcol/PVper closures but returns
#' the length-`np` gradient instead of the `np x np` R.  `dOiEst`/`tr28` are the
#' ESTIMATION-scale Omega-inverse derivatives (Cholesky params) -- `dOiEst[[k]] =
#' dOmega^{-1}/d(theta_omega_k)`, `tr28[k] = 0.5*tr(dOmega^{-1}_k %*% Omega)` --
#' NOT the variance-scale `omd` the cov uses for reported SEs.  Param order:
#' `nth` theta, `nsg` sigma, then `length(dOiEst)` Omega (Cholesky) params.
#' @noRd
.foceiAnalyticSubjectGrad <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, dOiEst, tr28,
                                      ndir = neta, dirTh = seq_len(nth), Oi = solve(Om)) {
  tr <- function(M) sum(diag(M))
  a <- E$a; A <- E$A; f <- E$f; y <- E$y
  evf <- function(e) ef$ev(e, f, y)
  rd <- list(r1 = evf(ef$sc$r1), r2 = evf(ef$sc$r2))
  pf <- list(p = evf(ef$sc$p), p1 = evf(ef$sc$p1))
  nom <- length(dOiEst)
  np <- nth + nsg + nom
  ei <- seq_len(neta); di <- seq_len(ndir)
  ae <- a[, ei, drop = FALSE]
  # exact inner Hessian H = d2l/deta2 (drives the EBE implicit derivative etaP);
  # Laplace determinant Hessian Ht = sum(p a a) + Omega^-1 (Eq 14).
  H <- Oi; for (l in ei) for (m in ei) H[l, m] <- H[l, m] + sum(rd$r2 * a[, l] * a[, m] + rd$r1 * A[, l, m])
  HiM <- solve(H)
  N <- matrix(0, neta, ndir); for (l in ei) for (d in di) N[l, d] <- sum(rd$r2 * a[, l] * a[, d] + rd$r1 * A[, l, d])
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(pf$p * a[, l] * a[, m]); Hti <- solve(Ht)
  ouAA <- function(v) { M <- matrix(0, neta, neta); for (l in ei) for (m in ei) M[l, m] <- sum(v * a[, l] * a[, m]); M }
  # dHt/d(direction s): eta-directions (s in ei) give the moving-mode dHt/d eta.
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(pf$p1 * a[, s] * a[, l] * a[, m] + pf$p * A[, l, s] * a[, m] + pf$p * a[, l] * A[, m, s]); D })
  typ <- function(p) if (p <= nth) "th" else if (p <= nth + nsg) "sg" else "om"
  sgi <- function(p) sgVar[p - nth]
  omc <- function(p) p - nth - nsg
  PVper <- function(p) lapply(ef$per[[sgi(p)]], evf)
  # M_p = d2l/(deta dp): th -> N[,dir]; sg -> a' (d2rho/df dsig); om -> dOmega^-1 eta
  Mcol <- function(p) { t <- typ(p)
    if (t == "th") return(N[, dirTh[p]]); if (t == "sg") return(as.numeric(crossprod(ae, PVper(p)$rf)))
    as.numeric(dOiEst[[omc(p)]] %*% ehat) }
  # explicit dHt/dp: th -> dHtD[dir]; sg -> outer(dp/dsig); om -> dOmega^-1
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") return(dHtD[[dirTh[p]]]); if (t == "sg") return(ouAA(PVper(p)$ps))
    dOiEst[[omc(p)]] }
  # explicit dPhi/dp (envelope): th -> sum(r1 df/dth); sg -> sum(drho/dsig);
  # om -> 0.5 eta' dOmega^-1 eta - 0.5 tr(dOmega^-1 Omega)
  dPhi_p <- function(p) { t <- typ(p)
    if (t == "th") return(sum(rd$r1 * a[, dirTh[p]]))
    if (t == "sg") return(sum(PVper(p)$rs))
    0.5 * as.numeric(t(ehat) %*% dOiEst[[omc(p)]] %*% ehat) - tr28[omc(p)] }
  # EBE parameter sensitivity etaP[,p] = d eta*/dp = -H^-1 M_p (Almquist Eq 46)
  etaP <- vapply(seq_len(np), function(p) as.numeric(-HiM %*% Mcol(p)), numeric(neta))
  if (neta == 1L) etaP <- matrix(etaP, nrow = 1L)
  g <- numeric(np)
  for (p in seq_len(np)) {
    # d(log|Ht|)/dp = tr(Ht^-1 (dHt/dp|explicit + sum_l dHt/deta_l * etaP[l,p]))
    dHtStar <- dHt_p(p)
    for (l in ei) dHtStar <- dHtStar + etaP[l, p] * dHtD[[l]]
    g[p] <- 2 * dPhi_p(p) + tr(Hti %*% dHtStar)      # d(OFV_i)/dp
  }
  list(g = g, etaP = etaP)                            # etaP = d eta*/d p (Almquist Eq 46/48)
}

#' C++/Armadillo port of `.foceiAnalyticSubjectGrad`: evaluates the per-observation
#' error coefficients in R (cheap, vectorized) and does the O(neta^2*nobs) tensor
#' contractions in `foceiSubjectGradFocei_`.  Same signature/return as the R oracle.
#' @noRd
.foceiAnalyticSubjectGradCpp <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, dOiEst, tr28,
                                         ndir = neta, dirTh = seq_len(nth), Oi = solve(Om)) {
  a <- E$a; A <- E$A; f <- E$f; y <- E$y
  nobs <- length(f)
  # A constant (f-independent) derivative -- e.g. d2rho/df2 = 1/sa^2 for a pure
  # additive error -- evaluates to a length-1 scalar; the C++ kernel indexes every
  # coefficient per observation, so recycle each to length nobs.
  .n <- function(e) rep_len(as.numeric(ef$ev(e, f, y)), nobs)
  r1 <- .n(ef$sc$r1); r2 <- .n(ef$sc$r2); pp <- .n(ef$sc$p); p1 <- .n(ef$sc$p1)
  perRf <- matrix(0, nobs, nsg); perPs <- matrix(0, nobs, nsg); perRs <- matrix(0, nobs, nsg)
  if (nsg > 0L) for (j in seq_len(nsg)) {
    perRf[, j] <- .n(ef$per[[sgVar[j]]]$rf)
    perPs[, j] <- .n(ef$per[[sgVar[j]]]$ps)
    perRs[, j] <- .n(ef$per[[sgVar[j]]]$rs)
  }
  nom <- length(dOiEst)
  dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
  foceiSubjectGradFocei_(a, A, r1, r2, pp, p1, perRf, perPs, perRs, as.numeric(ehat), Oi,
                         dOiCube, if (nom > 0L) as.numeric(tr28) else numeric(0),
                         neta, nth, nsg, nom, as.integer(dirTh))
}

#' Per-subject first-derivative (outer-gradient) contribution for FOCE
#' (interaction=0).  Unlike FOCEI, the EBE eta* stationarizes the interaction-free
#' inner problem S_FOCE = sum(q a) + Omega^-1 eta = 0 (NOT the full Laplace
#' objective), so the envelope shortcut fails: the gradient is the general total
#' derivative `dF/dp = dF/dp|explicit + F_eta . eta_p`, with `eta_p = -Hf^-1 S_p`
#' (Hf the FOCE inner Hessian) and `F_eta = Phi_eta + 0.5 d log|Ht|/d eta` nonzero.
#' `foceType=0` ("nonmem") freezes the variance R0 at the eta=0 population
#' prediction f0 (from E0), adding the a0-chain corrections; `foceType=1` ("foce+")
#' keeps the live conditional R.  Omega block uses the estimation-scale dOiEst/tr28.
#' @noRd
.foceiAnalyticSubjectGradFoce <- function(E, ehat, Om, ef, neta, nth, nsg, sgVar, dOiEst, tr28,
                                          ndir = neta, dirTh = seq_len(nth), Oi = solve(Om),
                                          E0 = NULL, foceType = 0L) {
  tr <- function(M) sum(diag(M))
  a <- E$a; A <- E$A; f <- E$f; y <- E$y
  .fc <- if (identical(as.integer(foceType), 1L)) ef$focePlus else ef$foce
  f0 <- if (!is.null(E0)) E0$f else f
  evf <- function(e) ef$ev(e, f, y, f0)
  rd <- list(r1 = evf(.fc$sc$r1))                                       # full rho (Phi), R0
  qd <- list(q0 = evf(.fc$sc$q0), q1 = evf(.fc$sc$q1))                  # FOCE inner gradient coef
  pf <- list(p = evf(.fc$sc$pF), p1 = evf(.fc$sc$pF1))                  # FOCE determinant p = 1/R0
  nom <- length(dOiEst); np <- nth + nsg + nom
  ei <- seq_len(neta); di <- seq_len(ndir); ae <- a[, ei, drop = FALSE]
  .dirOf <- function(p) dirTh[p]
  .cf0 <- isTRUE(.fc$dependsF0) && !is.null(E0)
  a0 <- if (.cf0) E0$a else NULL
  fq <- if (.cf0) list(qf0 = evf(.fc$f0$qf0), pFf0 = evf(.fc$f0$pFf0), rhof0 = evf(.fc$f0$rhof0)) else NULL
  # Phi_eta (nonzero at the FOCE eta*), FOCE inner Hessian Hf, its Nf, determinant Ht
  gPhi <- as.numeric(Oi %*% ehat); for (l in ei) gPhi[l] <- gPhi[l] + sum(rd$r1 * a[, l])
  Hf <- Oi; for (l in ei) for (m in ei) Hf[l, m] <- Hf[l, m] + sum(qd$q1 * a[, l] * a[, m] + qd$q0 * A[, l, m])
  HfInv <- solve(Hf)
  Nf <- matrix(0, neta, ndir); for (l in ei) for (d in di) Nf[l, d] <- sum(qd$q1 * a[, l] * a[, d] + qd$q0 * A[, l, d])
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(pf$p * a[, l] * a[, m]); Hti <- solve(Ht)
  ouAA <- function(v) { M <- matrix(0, neta, neta); for (l in ei) for (m in ei) M[l, m] <- sum(v * a[, l] * a[, m]); M }
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(pf$p1 * a[, s] * a[, l] * a[, m] + pf$p * A[, l, s] * a[, m] + pf$p * a[, l] * A[, m, s]); D })
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtD[[l]]), numeric(1))
  typ <- function(p) if (p <= nth) "th" else if (p <= nth + nsg) "sg" else "om"
  sgi <- function(p) sgVar[p - nth]; omc <- function(p) p - nth - nsg
  # S_p = dS_FOCE/dp (EBE stationarity), for eta_p = -Hf^-1 S_p
  McolEBE <- function(p) { t <- typ(p)
    if (t == "th") { v <- Nf[, .dirOf(p)]
      if (.cf0) v <- v + as.numeric(crossprod(ae, fq$qf0 * a0[, .dirOf(p)])); return(v) }
    if (t == "sg") return(as.numeric(crossprod(ae, evf(.fc$per[[sgi(p)]]$qs))))
    as.numeric(dOiEst[[omc(p)]] %*% ehat) }
  # explicit dHt/dp (p = 1/R0): th adds pFf0*a0 chain under foceType=0
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") { D <- dHtD[[.dirOf(p)]]
      if (.cf0) { d <- .dirOf(p); for (m in ei) for (n in ei) D[m, n] <- D[m, n] + sum(fq$pFf0 * a0[, d] * a[, m] * a[, n]) }
      return(D) }
    if (t == "sg") return(ouAA(evf(.fc$per[[sgi(p)]]$psF)))
    dOiEst[[omc(p)]] }
  # explicit dPhi/dp (fixed eta): th uses full r1 (+ rhof0*a0 chain); sg pure d(rho)/d(sig);
  # om the prior 0.5 eta' dOmega^-1 eta - 0.5 tr(dOmega^-1 Omega)
  dPhiExplicit <- function(p) { t <- typ(p)
    if (t == "th") { v <- sum(rd$r1 * a[, .dirOf(p)])
      if (.cf0) v <- v + sum(fq$rhof0 * a0[, .dirOf(p)]); return(v) }
    if (t == "sg") return(sum(evf(.fc$per[[sgi(p)]]$rs)))
    0.5 * as.numeric(t(ehat) %*% dOiEst[[omc(p)]] %*% ehat) - tr28[omc(p)] }
  etaP <- vapply(seq_len(np), function(p) as.numeric(-HfInv %*% McolEBE(p)), numeric(neta))
  if (neta == 1L) etaP <- matrix(etaP, nrow = 1L)
  g <- numeric(np)
  for (p in seq_len(np)) {
    # d(OFV_i)/dp = 2 (dPhi/dp|explicit + 0.5 tr(Ht^-1 dHt/dp|explicit) + (Phi_eta + Cen).eta_p)
    g[p] <- 2 * (dPhiExplicit(p) + 0.5 * tr(Hti %*% dHt_p(p)) + sum((gPhi + Cen) * etaP[, p]))
  }
  list(g = g, etaP = etaP)                            # etaP = d eta*/d p (Almquist Eq 46/48)
}

#' Solve the augmented model for ALL subjects in one rxode2 population solve
#' (internally C++ + OpenMP-threaded), instead of one R solve per subject -- the
#' per-subject solve loop otherwise dominates the gradient (~87%).  Per-subject
#' etas travel as an ID-keyed params data.frame.  Returns a per-subject list of
#' `list(f, a, A)` (same shape as `.foceiAnalyticSolveFA`), or `NULL` on failure.
#' FOCEI only (interaction=1, EBEs at the stored values -- no per-subject re-solve).
#' @noRd
.foceiAnalyticSolveAll <- function(am, thv, ebes, ids, data, obsTimes, tol = 1e-10) {
  dirs <- am$dirs; nd <- length(dirs); neta <- ncol(ebes)
  etav <- paste0("ETA_", seq_len(neta), "_")
  pars <- data.frame(ID = ids)
  for (k in seq_len(neta)) pars[[etav[k]]] <- ebes[, k]
  for (.nm in names(thv)) pars[[.nm]] <- thv[[.nm]]
  .sol <- tryCatch(withCallingHandlers(
    as.data.frame(rxode2::rxSolve(am$augMod, params = pars, events = data,
                                  returnType = "data.frame", atol = tol, rtol = tol)),
    warning = function(w) invokeRestart("muffleWarning")), error = function(e) NULL)
  if (is.null(.sol) || !all(c("rx_predf_", paste0("rx_f1_", dirs)) %in% names(.sol))) return(NULL)
  .idcol <- if ("id" %in% names(.sol)) .sol$id else .sol$ID
  .byIdSol <- split(seq_len(nrow(.sol)), as.character(.idcol))
  # rxode2 may label the solve id by the input ID or renumber 1..nsub; accept either.
  Es <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    .ri <- .byIdSol[[as.character(ids[i])]]
    if (is.null(.ri)) .ri <- .byIdSol[[as.character(i)]]
    if (is.null(.ri)) return(NULL)
    .di <- .sol[.ri, , drop = FALSE]
    .di <- .di[.di$time %in% obsTimes[[i]], , drop = FALSE]
    if (nrow(.di) != length(obsTimes[[i]])) return(NULL)
    no <- nrow(.di)
    a <- matrix(vapply(dirs, function(q) .di[[paste0("rx_f1_", q)]], numeric(no)), no, nd)
    A <- array(0, c(no, nd, nd))
    for (r in seq_len(nrow(am$P2))) {
      .ii <- match(am$P2$i[r], dirs); .jj <- match(am$P2$j[r], dirs)
      .v2 <- .di[[paste0("rx_f2_", am$P2$i[r], "_", am$P2$j[r])]]
      A[, .ii, .jj] <- .v2; A[, .jj, .ii] <- .v2
    }
    Es[[i]] <- list(f = .di$rx_predf_, a = a, A = A)
    if (isTRUE(am$hasRvar)) Es[[i]] <- c(Es[[i]], .foceiReadRvar(.di, am, no))
  }
  Es
}

#' Shared core: analytic natural-scale outer gradient of the FOCEI objective
#' (OFV = -2*logLik) over structural theta + residual sigma + Omega (Cholesky
#' scale), summed over subjects.  Gathering (theta/eta/data/omega) is done by the
#' live-env and post-fit callers; this runs the augmented solve + assembly.
#' Returns a NAMED numeric vector or `NULL`.
#' @noRd
.foceiAnalyticGradCore <- function(ui, th, ebes, ids, data, Om, ef, .dir, dOiEst, tr28,
                                   omNames, solveTol, interaction = 1L, foceType = 0L,
                                   startedEnv = NULL, am = NULL) {
  neta <- ncol(ebes); Oi <- solve(Om)
  thStruct <- .dir$thStruct; dirs <- .dir$dirs; dirTh <- .dir$dirTh; ndir <- .dir$ndir; nth <- .dir$nth
  nsg <- length(ef$sgVar); nom <- length(dOiEst); np <- nth + nsg + nom
  etav <- paste0("ETA_", seq_len(neta), "_")
  .foce <- identical(as.integer(interaction), 0L)
  # The augmented model depends only on the model + direction set (fixed for a
  # fit), NOT on theta/eta/omega; the symbolic .rxSens build dominates each
  # gradient (~63%), so the live path passes a cached `am` (built once per fit).
  if (is.null(am)) am <- .foceiAnalyticAugModelDirs(ui, dirs)
  if (is.null(am) || am$ndir != ndir) return(NULL)
  .byId <- split(data, as.character(data$ID))
  .idCode <- if (is.factor(ids)) as.integer(ids) else match(ids, sort(unique(ids)))
  if (!is.null(startedEnv)) assign(".analyticStarted", TRUE, startedEnv)
  .obsT <- lapply(seq_along(ids), function(i) {
    .s <- .byId[[as.character(.idCode[i])]]; .s$TIME[.s$EVID == 0]
  })
  # FOCE: the per-subject EBE re-solve (Newton) + optional eta=0 population solve are
  # inherently per-subject; collect the re-solved etas (and E0) here.  FOCEI EBEs are
  # fixed.  The FINAL a/A solve is then batched for ALL subjects in one rxode2
  # population solve (C++ + OpenMP), which otherwise dominates the gradient cost.
  etaSolve <- ebes
  E0List <- vector("list", length(ids))
  if (.foce) {
    for (i in seq_along(ids)) {
      s <- .byId[[as.character(.idCode[i])]]; if (is.null(s) || nrow(s) == 0L) return(NULL)
      obs <- s[s$EVID == 0, , drop = FALSE]
      E0 <- NULL
      if (identical(as.integer(foceType), 0L) && isTRUE(ef$foce$dependsF0)) {
        E0 <- .foceiAnalyticSolveFA(am, c(th, setNames(rep(0, neta), etav)), s, obs$TIME, tol = solveTol)
        if (is.null(E0)) return(NULL)
      }
      E0List[[i]] <- E0
      eta0 <- .foceiAnalyticFoceEbe(am, th, ebes[i, ], s, obs$TIME, obs$DV, etav, ef, Oi, neta, solveTol,
                                    f0 = if (is.null(E0)) NULL else E0$f, foceType = foceType)
      if (is.null(eta0)) return(NULL)
      etaSolve[i, ] <- eta0
    }
  }
  .EsAll <- .foceiAnalyticSolveAll(am, th, etaSolve, .idCode, data, .obsT, solveTol)
  if (is.null(.EsAll)) return(NULL)
  g <- numeric(np)
  etaPList <- vector("list", length(ids))            # per-subject d eta*/d p (Eq 48)
  for (i in seq_along(ids)) {
    obs <- .byId[[as.character(.idCode[i])]]
    obs <- obs[obs$EVID == 0, , drop = FALSE]
    E <- .EsAll[[i]]; if (is.null(E)) return(NULL)
    if (isTRUE(ef$canVanish)) {
      .fa <- abs(E$f); if (any(!is.finite(.fa)) || min(.fa) < 1e-6 * max(.fa)) return(NULL)
    }
    E$y <- obs$DV
    gi <- tryCatch(
      if (.foce) .foceiAnalyticSubjectGradFoce(E, etaSolve[i, ], Om, ef, neta, nth, nsg, ef$sgVar,
                                               dOiEst, tr28, ndir = ndir, dirTh = dirTh, Oi = Oi,
                                               E0 = E0List[[i]], foceType = foceType)
      else .foceiAnalyticSubjectGradCpp(E, etaSolve[i, ], Om, ef, neta, nth, nsg, ef$sgVar,
                                        dOiEst, tr28, ndir = ndir, dirTh = dirTh, Oi = Oi),
      error = function(e) NULL)
    if (is.null(gi) || !all(is.finite(gi$g)) || !all(is.finite(gi$etaP))) return(NULL)
    g <- g + gi$g
    etaPList[[i]] <- gi$etaP
  }
  names(g) <- c(thStruct, ef$sgName, omNames)
  list(g = g, etaP = etaPList, ids = ids)
}

#' Common scope gates + error/direction/omega setup shared by the live and
#' post-fit gradient paths.  Returns a list of the assembled pieces, or `NULL`
#' (out of scope).  `thVals` is the named converged theta vector.
#' @noRd
.foceiAnalyticGradSetup <- function(ui, thVals, Om) {
  if (!isTRUE(rxode2::rxGetControl(ui, "fast", FALSE))) return(NULL)
  if (!.hasRxSens()) return(NULL)
  # bounded parameter transforms are corrected on a different (natural) scale
  if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L) return(NULL)
  # tad/podo/tafd/tlast/tfirst/dosenum are functions of time and the dose record
  # only (no eta/theta dependence), so rxode2 treats them as zero-derivative
  # constants in the sensitivity expansion (.rxToSEDualVarFunction) -- they no
  # longer need to force the finite-difference fallback.
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(NULL)
  interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))            # 1 FOCEI / 0 FOCE
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(ui, "foceType", 0L)) else 0L
  # foce+ (foceType=1, live conditional R): the inner optimizer uses the
  # interaction-free gradient (eps*a/R) while the objective/likelihood uses the
  # LIVE R (which depends on eta), so its function and gradient are inconsistent
  # and the EBE does not stationarize S_FOCE=0 (observed max|S_FOCE|~1 even at a
  # converged fit).  The analytic gradient assumes a clean S_FOCE=0 stationary
  # point (as the covariance does), so it cannot match the finite-difference
  # gradient of that ill-defined EBE -- keep foce+ on the FD gradient.
  if (interaction == 0L && foceType == 1L) return(NULL)
  if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L) return(NULL)
  ef <- .foceiAnalyticErrFull(ui); if (is.null(ef)) return(NULL)
  ini <- ui$iniDf
  .map <- .foceiEtaThetaMap(ui)
  etaNames <- .map$etaNames; neta <- length(etaNames)
  if (neta == 0L) return(NULL)
  thetaForEta <- .map$thetaForEta
  if (length(.uiIovEnv$iovVars) > 0L) return(NULL)                                  # IOV out of Phase-1 scope
  .valc <- setNames(as.numeric(thVals[ef$sgName]), ef$sgVar)
  ef$ev <- local({ v <- .valc; function(expr, f, y, f0 = f) eval(expr, c(list(f = f, y = y, f0 = f0), as.list(v))) })
  if (any(.iniIsFixed(ini, thetaForEta))) return(NULL)
  keep <- !.iniIsFixed(ini, ef$sgName); ef$sgVar <- ef$sgVar[keep]; ef$sgName <- ef$sgName[keep]
  .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta)
  if (is.null(.dir)) return(NULL)
  .oe <- .foceiEstOmegaDeriv(ui, Om); if (is.null(.oe)) return(NULL)
  list(ef = ef, dir = .dir, dOiEst = .oe$dOi, tr28 = .oe$tr28, omNames = .oe$names,
       neta = neta, etaNames = etaNames, interaction = interaction, foceType = foceType)
}

#' Post-fit analytic natural-scale gradient for a fitted object (validation /
#' standalone).  Mirrors `.foceiCovAnalyticCalc`'s gathering.
#' @noRd
.foceiGradAnalyticCalc <- function(fit) {
  tryCatch({
    ui <- fit$finalUi
    Om <- fit$omega
    ini <- ui$iniDf
    thRows <- ini[!is.na(ini$ntheta), , drop = FALSE]
    thRows <- thRows[order(thRows$ntheta), , drop = FALSE]
    .thv <- fit$theta[thRows$name]; if (anyNA(.thv)) .thv <- thRows$est
    thVals <- setNames(as.numeric(.thv), thRows$name)
    st <- .foceiAnalyticGradSetup(ui, thVals, Om)
    if (is.null(st)) return(NULL)
    th <- setNames(as.numeric(thVals), paste0("THETA_", seq_along(thVals), "_"))
    ebes <- as.matrix(fit$eta[, st$etaNames, drop = FALSE])
    if (!is.null(fit$dataSav$CENS) && any(fit$dataSav$CENS != 0, na.rm = TRUE)) return(NULL)
    .r <- .foceiAnalyticGradCore(ui, th, ebes, fit$eta$ID, fit$dataSav, Om, st$ef, st$dir,
                                 st$dOiEst, st$tr28, st$omNames, .foceiAnalyticSolveTol(ui),
                                 interaction = st$interaction, foceType = st$foceType)
    if (is.null(.r)) return(NULL)
    .r$g                                              # named natural-scale gradient (validation/tests)
  }, error = function(e) NULL)
}

#' Live-fit analytic natural-scale gradient from the C++ callback env `e`.
#' @noRd
.foceiAnalyticGradFocei <- function(e) {
  tryCatch({
    ui <- get("ui", e)
    Om <- get("omega", e)
    thNames <- get("thetaNames", e)
    thVals  <- get("theta", e)$theta; names(thVals) <- thNames
    st <- .foceiAnalyticGradSetup(ui, thVals, Om)
    if (is.null(st)) return(NULL)
    th <- setNames(as.numeric(thVals[thNames]), paste0("THETA_", seq_along(thNames), "_"))
    etaObf <- get("etaObf", e)
    ebes <- as.matrix(etaObf[, paste0("ETA[", seq_len(st$neta), "]"), drop = FALSE])
    data <- get("dataSav", e)
    if (!is.null(data$CENS) && any(data$CENS != 0, na.rm = TRUE)) return(NULL)
    # The augmented model is the persistent `..outer` sibling of the inner model.
    # Prefer the copy built at model-setup time and qs2-cached in foceiModel$outer
    # (reconstruct am from the top-level compiled model + outerMeta); fall back to
    # building it via rxUiGet.foceiOuter.  Cached on the fit env either way.
    am <- if (exists(".foceiGradAug", e, inherits = FALSE)) get(".foceiGradAug", e) else NULL
    if (is.null(am)) {
      .fm <- tryCatch(ui$foceiModel, error = function(e) NULL)
      if (!is.null(.fm) && inherits(.fm$outer, "rxode2") && !is.null(.fm$outerMeta)) {
        am <- c(list(augMod = .fm$outer), .fm$outerMeta)
      } else {
        am <- ui$foceiOuter
      }
      if (!is.null(am) && inherits(am$augMod, "rxode2")) assign(".foceiGradAug", am, envir = e)
    }
    .foceiAnalyticGradCore(ui, th, ebes, etaObf$ID, data, Om, st$ef, st$dir,
                           st$dOiEst, st$tr28, st$omNames, .foceiAnalyticSolveTol(ui),
                           interaction = st$interaction, foceType = st$foceType,
                           startedEnv = e, am = am)
  }, error = function(e) NULL)
}

#' Direction set for the augmented outer-gradient model, computed from the UI
#' alone (does not depend on theta/eta values): one direction per eta plus one per
#' non-mu-referenced structural theta.  `NULL` if out of analytic scope.
#' @noRd
.foceiOuterDirs <- function(ui) {
  if (!.hasRxSens()) return(NULL)
  if (!is.null(ui$boundedTransforms) && length(ui$boundedTransforms) > 0L) return(NULL)
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(NULL)
  ef <- .foceiAnalyticErrFull(ui); if (is.null(ef)) return(NULL)
  .map <- .foceiEtaThetaMap(ui); neta <- length(.map$etaNames)
  if (neta == 0L) return(NULL)
  if (length(.uiIovEnv$iovVars) > 0L) return(NULL)
  .foceiAnalyticDirections(ui$iniDf, .map$thetaForEta, ef$sgName, neta)
}

#' Build the augmented outer-gradient sensitivity model (compiled model + `dirs` +
#' `P2`) for a UI.  This is the persistent `..outer` sibling of the inner model:
#' it depends only on the model + direction set (NOT theta/eta/omega), so it is
#' built once during model setup (via `rxUiGet.foceiModel`/`foceModel`, which
#' qs2-cache the whole model list) and reused across every outer-gradient call.
#' Callable independently as `ui$foceiOuter`.  `NULL` when out of analytic scope
#' (the gradient then falls back to finite differences).
#' @export
rxUiGet.foceiOuter <- function(x, ...) {
  .ui <- x[[1]]
  if (!isTRUE(rxode2::rxGetControl(.ui, "fast", FALSE))) return(NULL)
  interaction <- as.integer(rxode2::rxGetControl(.ui, "interaction", 1L))
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(.ui, "foceType", 0L)) else 0L
  if (interaction == 0L && foceType == 1L) return(NULL)              # foce+ stays on FD
  if (as.integer(rxode2::rxGetControl(.ui, "nAGQ", 1L)) > 1L) return(NULL)
  .dir <- .foceiOuterDirs(.ui); if (is.null(.dir)) return(NULL)
  .foceiAnalyticAugModelDirs(.ui, .dir$dirs)
}
attr(rxUiGet.foceiOuter, "rstudio") <- emptyenv()

#' Estimation-scale (Cholesky) Omega-inverse derivatives for the outer gradient's
#' Omega block, from rxode2's `rxSymInvCholEnvCalculate` (rxSymInv.R d.omegaInv /
#' tr.28).  Returns `list(dOi = list(dOmega^{-1}/d theta_omega_k), tr28 = 0.5 *
#' tr(dOmega^{-1}_k Omega), names = <omega parameter names>)`, or `NULL`.
#' @noRd
.foceiEstOmegaDeriv <- function(ui, Om) {
  tryCatch({
    # Build the rxSymInvChol env from the current Omega with the SAME diagonal
    # transform the optimizer uses, so the Cholesky parameter order/scale matches
    # op_focei's Omega slots (rxUiGet.focei builds env$rxInv the same way).
    .diagXform <- rxode2::rxGetControl(ui, "diagXform", "sqrt")
    .rxInv <- rxode2::rxSymInvCholCreate(mat = Om, diag.xform = .diagXform)
    # `$.rxSymInvCholEnv` dispatches to the C rxSymInvCholEnvCalculate; d.omegaInv
    # is the list of dOmega^-1/d(chol theta_k), tr.28 = 0.5*tr(dOmega^-1_k Omega).
    .dOi <- .rxInv$d.omegaInv
    .tr28 <- .rxInv$tr.28
    if (is.null(.dOi) || is.null(.tr28)) return(NULL)
    list(dOi = .dOi, tr28 = as.numeric(.tr28),
         names = paste0("om.chol.", seq_along(.dOi)))
  }, error = function(e) NULL)
}

#' Analytic outer gradient of the FOCE/FOCEI objective (Almquist 2015 Eq 23).
#'
#' Called from C++ (`analyticOuterGrad`) when `foceiControl(fast=TRUE)`.  The fit
#' environment `e` carries the current theta/eta state; the C++ caller has already
#' evaluated the objective at the current theta so the inner solutions (eta*) are
#' current.
#'
#' @param e focei fit environment
#' @return numeric gradient vector (length npars, model-theta scale), or `NULL`
#'   to fall back to the finite-difference gradient
#' @noRd
.foceiCalcGradAnalytic <- function(e) {
  tryCatch({
    .r <- .foceiAnalyticGradFocei(e)
    if (is.null(.r) || is.null(.r$g) || !all(is.finite(.r$g))) return(NULL)
    g <- .r$g; gn <- names(g)
    # op_focei parameter order is fullTheta = [non-fixed thetas by ntheta order |
    # omega Cholesky params] (inner.cpp fullTheta layout).  Map the named
    # natural-scale gradient onto that order; the C++ hook length-checks against
    # op_focei.npars and falls back to FD on any mismatch.
    thNames <- get("thetaNames", e)
    thOrder <- thNames[thNames %in% gn]                 # non-fixed structural+sigma thetas, in order
    omIdx <- grep("^om\\.chol\\.", gn)
    parOrder <- c(match(thOrder, gn), omIdx)            # gradient/etaP column order -> npars order
    gvec <- g[parOrder]
    if (anyNA(gvec) || !all(is.finite(gvec))) return(NULL)
    # Stash the per-subject d eta*/d p (etaP), columns reordered to the npars order,
    # for the C++ Eq-48 extrapolation (mceta=-2/-1).  neta x npars x nsub array,
    # aligned to etaObf$ID order (same as the fit's per-subject solve order).
    .etaP <- .r$etaP
    if (length(.etaP) > 0L && !any(vapply(.etaP, is.null, logical(1)))) {
      .neta <- nrow(.etaP[[1]]); .np <- length(parOrder); .ns <- length(.etaP)
      .arr <- array(NA_real_, c(.neta, .np, .ns))
      for (i in seq_len(.ns)) .arr[, , i] <- .etaP[[i]][, parOrder, drop = FALSE]
      if (all(is.finite(.arr))) assign(".foceiGradEtaP", .arr, envir = e)
      else if (exists(".foceiGradEtaP", e, inherits = FALSE)) rm(".foceiGradEtaP", envir = e)
    }
    as.numeric(gvec)
  }, error = function(e) NULL)
}
