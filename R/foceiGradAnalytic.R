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

#' (f,R) FOCEI per-subject outer gradient: the general form that reads the residual
#' variance R and its sensitivities (aR/AR) from the solve, treating the transformed
#' prediction f and the variance R as INDEPENDENT solved quantities.  The rho(f,R,y)
#' partials are model-independent closed forms, so ANY variance structure works.  A
#' residual sigma is a pseudo-direction with df/dsigma=0 and dR/dsigma=E$Rsig; the
#' omega block uses the estimation-scale dOiEst/tr28.  `dirTh` maps each structural
#' theta param to its direction column; `sigCol` maps each sigma param to its E$Rsig
#' column.  Reduces exactly to `.foceiAnalyticSubjectGrad` when R=R(f).
#' @noRd
.foceiAnalyticSubjectGradFR <- function(E, ehat, Om, neta, nth, nsg, dirTh, sigCol, dOiEst, tr28,
                                        ndir, Oi = solve(Om)) {
  tr <- function(M) sum(diag(M))
  f <- E$f; y <- E$y; R <- E$R; a <- E$a; A <- E$A; aR <- E$aR; AR <- E$AR
  Rsig <- E$Rsig; RsigDir <- E$RsigDir
  res <- y - f
  rf <- -res / R; rR <- 0.5 * (1 / R - res^2 / R^2)                     # rho_f, rho_R
  rff <- 1 / R; rfR <- res / R^2; rRR <- 0.5 * (-1 / R^2 + 2 * res^2 / R^3)  # rho_ff/_fR/_RR
  eff <- 1 / R; eRR <- 0.5 / R^2                                        # E[rho_ff], E[rho_RR] at y=f
  nom <- length(dOiEst); np <- nth + nsg + nom
  ei <- seq_len(neta); di <- seq_len(ndir)
  # exact inner Hessian H (eta x eta) and N (eta x direction) via the (f,R) contraction
  H <- Oi; N <- matrix(0, neta, ndir)
  for (l in ei) {
    for (m in ei) H[l, m] <- H[l, m] + sum(rff * a[, l] * a[, m] + rfR * (a[, l] * aR[, m] + aR[, l] * a[, m]) +
                                             rRR * aR[, l] * aR[, m] + rf * A[, l, m] + rR * AR[, l, m])
    for (d in di) N[l, d] <- sum(rff * a[, l] * a[, d] + rfR * (a[, l] * aR[, d] + aR[, l] * a[, d]) +
                                   rRR * aR[, l] * aR[, d] + rf * A[, l, d] + rR * AR[, l, d])
  }
  HiM <- solve(H)
  # eta x sigma block (df/dsigma=0): Nsg[l,k]
  Nsg <- matrix(0, neta, nsg)
  if (nsg > 0L) for (l in ei) for (k in seq_len(nsg)) { .c <- sigCol[k]
    Nsg[l, k] <- sum(rfR * a[, l] * Rsig[, .c] + rRR * aR[, l] * Rsig[, .c] + rR * RsigDir[, l, .c]) }
  # Laplace determinant Hessian Ht = Oi + sum(E[rho_ff] a a + E[rho_RR] aR aR)
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(eff * a[, l] * a[, m] + eRR * aR[, l] * aR[, m])
  Hti <- solve(Ht)
  # dHt/d(direction s) (moving mode s=eta AND explicit theta columns): de_ff/ds=-aR_s/R^2
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(-aR[, s] / R^2 * a[, l] * a[, m] + eff * (A[, l, s] * a[, m] + a[, l] * A[, m, s]) +
                     -aR[, s] / R^3 * aR[, l] * aR[, m] + eRR * (AR[, l, s] * aR[, m] + aR[, l] * AR[, m, s])); D })
  # dHt/dsigma_k (df/dsigma=0 so no d2f term)
  dHtSg <- if (nsg > 0L) lapply(seq_len(nsg), function(k) { .c <- sigCol[k]; D <- matrix(0, neta, neta)
    for (l in ei) for (m in ei)
      D[l, m] <- sum(-Rsig[, .c] / R^2 * a[, l] * a[, m] - Rsig[, .c] / R^3 * aR[, l] * aR[, m] +
                       eRR * (RsigDir[, l, .c] * aR[, m] + aR[, l] * RsigDir[, m, .c])); D }) else list()
  typ <- function(p) if (p <= nth) "th" else if (p <= nth + nsg) "sg" else "om"
  omc <- function(p) p - nth - nsg
  Mcol <- function(p) { t <- typ(p)
    if (t == "th") return(N[, dirTh[p]]); if (t == "sg") return(Nsg[, p - nth])
    as.numeric(dOiEst[[omc(p)]] %*% ehat) }
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") return(dHtD[[dirTh[p]]]); if (t == "sg") return(dHtSg[[p - nth]])
    dOiEst[[omc(p)]] }
  dPhi_p <- function(p) { t <- typ(p)
    if (t == "th") return(sum(rf * a[, dirTh[p]] + rR * aR[, dirTh[p]]))
    if (t == "sg") return(sum(rR * Rsig[, sigCol[p - nth]]))
    0.5 * as.numeric(t(ehat) %*% dOiEst[[omc(p)]] %*% ehat) - tr28[omc(p)] }
  etaP <- vapply(seq_len(np), function(p) as.numeric(-HiM %*% Mcol(p)), numeric(neta))
  if (neta == 1L) etaP <- matrix(etaP, nrow = 1L)
  g <- numeric(np)
  for (p in seq_len(np)) {
    dHtStar <- dHt_p(p); for (l in ei) dHtStar <- dHtStar + etaP[l, p] * dHtD[[l]]
    g[p] <- 2 * dPhi_p(p) + tr(Hti %*% dHtStar)
  }
  list(g = g, etaP = etaP)
}

#' (f,R) FOCE per-subject outer gradient: the general form for any variance structure.
#' FOCE evaluates the focep error structure with the variance frozen -- for `foceType=0`
#' ("nonmem") R0 and its sensitivities come from the eta=0 population solve `E0` with the
#' eta-direction sensitivities zeroed (frozen w.r.t. eta); for `foceType=1` ("foce+") R0
#' is the live conditional variance from the eta-hat solve `E`.  The inner problem is
#' interaction-free (S_FOCE = sum(q0 a) + Omega^-1 eta, q0 = -(y-f)/R0), so the envelope
#' shortcut fails and the gradient carries the moving mode via etaP = -Hf^-1 dS/dp.
#' @noRd
.foceiAnalyticSubjectGradFoceFR <- function(E, ehat, Om, neta, nth, nsg, dirTh, sigCol, dOiEst, tr28,
                                            ndir, Oi = solve(Om), E0 = NULL, foceType = 0L) {
  tr <- function(M) sum(diag(M))
  f <- E$f; y <- E$y; a <- E$a; A <- E$A
  ei <- seq_len(neta); di <- seq_len(ndir)
  .fp <- identical(as.integer(foceType), 1L) || is.null(E0)
  # variance source R0.  aRe drives the ETA-BLOCK (gPhi/dHtD/Ht): for foce+ it is the
  # live dR/ddir; for nonmem it is 0 (R0 frozen w.r.t. eta -- and w.r.t. the eta-hat f,
  # since R0 uses the separate eta=0 population prediction).  aRc drives the PARAMETER
  # COLUMNS (dR0/dtheta the theta-chain): live dR/ddir for foce+, else the eta=0 solve's
  # dR0/ddir -- NOT zeroed, so a mu-referenced theta (whose direction is an eta direction)
  # still gets its frozen-R0 theta sensitivity.
  if (.fp) { R0 <- E$R; aRe <- E$aR; aRc <- E$aR; R0sig <- E$Rsig }
  else { R0 <- E0$R; aRe <- matrix(0, nrow(a), ndir); aRc <- E0$aR; R0sig <- E0$Rsig }
  res <- y - f
  rho_f <- -res / R0; rho_R <- 0.5 * (1 / R0 - res^2 / R0^2)
  q0 <- -res / R0; q1 <- 1 / R0                                  # interaction-free inner gradient/curv
  nom <- length(dOiEst); np <- nth + nsg + nom
  ae <- a[, ei, drop = FALSE]
  # Phi_eta (nonzero at the FOCE eta*): rho_f df/deta + rho_R dR0/deta (aRe=0 for nonmem)
  gPhi <- as.numeric(Oi %*% ehat)
  for (l in ei) gPhi[l] <- gPhi[l] + sum(rho_f * a[, l] + rho_R * aRe[, l])
  # FOCE inner Hessian Hf (interaction-free), its Nf, and the determinant Ht = Oi + sum(a a/R0)
  Hf <- Oi; Nf <- matrix(0, neta, ndir); Ht <- Oi
  for (l in ei) {
    for (m in ei) { Hf[l, m] <- Hf[l, m] + sum(q1 * a[, l] * a[, m] + q0 * A[, l, m])
      Ht[l, m] <- Ht[l, m] + sum(a[, l] * a[, m] / R0) }
    for (d in di) Nf[l, d] <- sum(q1 * a[, l] * a[, d] + q0 * A[, l, d])
  }
  HfInv <- solve(Hf); Hti <- solve(Ht)
  # dHt/d(direction s) eta-block moving mode: d(1/R0)/ds a a + 1/R0 (A a + a A); aRe=0 (nonmem)
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(-aRe[, s] / R0^2 * a[, l] * a[, m] + (A[, l, s] * a[, m] + a[, l] * A[, m, s]) / R0); D })
  Cen <- vapply(ei, function(l) 0.5 * tr(Hti %*% dHtD[[l]]), numeric(1))
  typ <- function(p) if (p <= nth) "th" else if (p <= nth + nsg) "sg" else "om"
  omc <- function(p) p - nth - nsg
  # explicit dHt/dtheta a0-chain for nonmem (aRe=0 -> not already in dHtD): -aRc/R0^2 a a
  ouRc <- function(v) { M <- matrix(0, neta, neta); for (l in ei) for (m in ei) M[l, m] <- sum(v * a[, l] * a[, m]); M }
  # dS_FOCE/dp for etaP = -Hf^-1 S_p: th adds the R0-theta chain (res/R0^2) dR0/dtheta a
  McolEBE <- function(p) { t <- typ(p)
    if (t == "th") { d <- dirTh[p]; return(Nf[, d] + as.numeric(crossprod(ae, (res / R0^2) * aRc[, d]))) }
    if (t == "sg") return(as.numeric(crossprod(ae, (res / R0^2) * R0sig[, sigCol[p - nth]])))
    as.numeric(dOiEst[[omc(p)]] %*% ehat) }
  # explicit dHt/dp: th -> dHtD[dir] (+ a0-chain for nonmem); sg -> -R0sig/R0^2 a a; om -> dOmega^-1
  dHt_p <- function(p) { t <- typ(p)
    if (t == "th") { d <- dirTh[p]; D <- dHtD[[d]]
      if (!.fp) D <- D + ouRc(-aRc[, d] / R0^2); return(D) }
    if (t == "sg") return(ouRc(-R0sig[, sigCol[p - nth]] / R0^2))
    dOiEst[[omc(p)]] }
  # explicit dPhi/dp (fixed eta): rho_f df/dp + rho_R dR0/dp (columns use aRc)
  dPhiExplicit <- function(p) { t <- typ(p)
    if (t == "th") { d <- dirTh[p]; return(sum(rho_f * a[, d] + rho_R * aRc[, d])) }
    if (t == "sg") return(sum(rho_R * R0sig[, sigCol[p - nth]]))
    0.5 * as.numeric(t(ehat) %*% dOiEst[[omc(p)]] %*% ehat) - tr28[omc(p)] }
  etaP <- vapply(seq_len(np), function(p) as.numeric(-HfInv %*% McolEBE(p)), numeric(neta))
  if (neta == 1L) etaP <- matrix(etaP, nrow = 1L)
  g <- numeric(np)
  for (p in seq_len(np))
    g[p] <- 2 * (dPhiExplicit(p) + 0.5 * tr(Hti %*% dHt_p(p)) + sum((gPhi + Cen) * etaP[, p]))
  list(g = g, etaP = etaP)
}

#' C++/Armadillo port of `.foceiAnalyticSubjectGradFR`: does the O(neta^2*nobs) (f,R)
#' tensor contractions in `foceiSubjectGradFR_` (the rho(f,R,y) partials are computed in
#' C++ from f/y/R).  Same signature/return as the R FR assembler.
#' @noRd
.foceiAnalyticSubjectGradFRCpp <- function(E, ehat, Om, neta, nth, nsg, dirTh, sigCol, dOiEst, tr28,
                                           ndir, Oi = solve(Om), dvSens = matrix(0, length(E$f), 0L)) {
  nom <- length(dOiEst); nobs <- length(E$f)
  dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
  Rsig <- if (is.null(E$Rsig)) matrix(0, nobs, 0L) else E$Rsig
  RsigDir <- if (is.null(E$RsigDir)) array(0, c(nobs, ndir, 0L)) else E$RsigDir
  foceiSubjectGradFR_(E$a, E$A, E$aR, E$AR, Rsig, RsigDir, dvSens, E$f, E$y, E$R,
                      as.numeric(ehat), Oi, dOiCube, if (nom > 0L) as.numeric(tr28) else numeric(0),
                      neta, nth, nsg, nom, as.integer(dirTh), as.integer(sigCol))
}

#' C++/Armadillo port of `.foceiAnalyticSubjectGradFoceFR` (FOCE (f,R) outer gradient).
#' Resolves the frozen-R0 sensitivities (eta-block aRe: 0 for nonmem, live E$aR for foce+;
#' parameter columns aRc: E0$aR for nonmem, live E$aR for foce+; R0sig likewise) and calls
#' the kernel.  Matches `.foceiAnalyticSubjectGradFoceFR` exactly.
#' @noRd
.foceiAnalyticSubjectGradFoceFRCpp <- function(E, ehat, Om, neta, nth, nsg, dirTh, sigCol, dOiEst, tr28,
                                               ndir, Oi = solve(Om), E0 = NULL, foceType = 0L,
                                               dvSens = matrix(0, length(E$f), 0L)) {
  nobs <- length(E$f); nom <- length(dOiEst)
  .fp <- identical(as.integer(foceType), 1L) || is.null(E0)
  if (.fp) { R0 <- E$R; aRe <- E$aR; aRc <- E$aR; R0sig <- E$Rsig }
  else { R0 <- E0$R; aRe <- matrix(0, nobs, ndir); aRc <- E0$aR; R0sig <- E0$Rsig }
  if (is.null(R0sig)) R0sig <- matrix(0, nobs, nsg)
  dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
  foceiSubjectGradFoceFR_(E$a, E$A, aRe, aRc, R0sig, dvSens, E$f, E$y, R0, as.numeric(ehat), Oi,
                          dOiCube, if (nom > 0L) as.numeric(tr28) else numeric(0),
                          neta, nth, nsg, nom, as.integer(dirTh), as.integer(sigCol), as.integer(.fp))
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
    if (isTRUE(am$hasTrans))                        # both-sides transform: DV -> tbs(DV) scale
      Es[[i]]$trans <- list(yj = .di$rx_tyj_, lambda = .di$rx_tlambda_, low = .di$rx_tlow_, hi = .di$rx_thi_)
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
  lamDir <- .dir$lamDir; lamNames <- .dir$lamNames; lamIdx <- .dir$lamIdx
  nom <- length(dOiEst)
  etav <- paste0("ETA_", seq_len(neta), "_")
  .foce <- identical(as.integer(interaction), 0L)
  # The augmented model depends only on the model + direction set (fixed for a
  # fit), NOT on theta/eta/omega; the symbolic .rxSens build dominates each
  # gradient (~63%), so the live path passes a cached `am` (built once per fit).
  if (is.null(am)) am <- .foceiAnalyticAugModelDirs(ui, dirs)
  if (is.null(am) || am$ndir != ndir) return(NULL)
  # FOCEI (interaction=1) uses the general (f,R) assembly: the sigma set is EVERY
  # error parameter the variance depends on (am$sigTh), which covers any variance
  # Both FOCEI and FOCE assemble over am$sigTh (every error param the variance depends on)
  # via the (f,R) path, so any variance structure works for both.
  .sigTh <- am$sigTh
  .sgNameFR <- if (length(.sigTh)) ui$iniDf$name[match(.sigTh, ui$iniDf$ntheta)] else character(0)
  nsg <- length(.sigTh); sgNames <- .sgNameFR
  np <- nth + nsg + nom
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
      if (identical(as.integer(foceType), 0L) && isTRUE(ef$dependsF0)) {
        E0 <- .foceiAnalyticSolveFA(am, c(th, setNames(rep(0, neta), etav)), s, obs$TIME, tol = solveTol)
        if (is.null(E0)) return(NULL)
      }
      E0List[[i]] <- E0
      eta0 <- .foceiAnalyticFoceEbe(am, th, ebes[i, ], s, obs$TIME, obs$DV, etav,
                                    if (is.null(E0)) NULL else E0$R, Oi, neta, solveTol, foceType = foceType)
      if (is.null(eta0)) return(NULL)
      etaSolve[i, ] <- eta0
    }
  }
  .EsAll <- .foceiAnalyticSolveAll(am, th, etaSolve, .idCode, data, .obsT, solveTol)
  if (is.null(.EsAll)) return(NULL)
  g <- numeric(np)
  etaPList <- vector("list", length(ids))            # per-subject d eta*/d p (Eq 48)
  if (.foce) {
    # FOCE: assemble ALL subjects in ONE OpenMP C++ call (foceiGradAllFoceFR_).  The frozen
    # R0 sensitivities are resolved per subject in R (nonmem: aRe=0, aRc/R0/R0sig from the
    # eta=0 solve E0; foce+ / additive nonmem: all from the eta-hat solve E) then batched.
    nsub <- length(ids); nobsAll <- integer(nsub)
    .fpG <- identical(as.integer(foceType), 1L) || is.null(E0List[[1L]])   # uniform across subjects
    for (i in seq_len(nsub)) {
      E <- .EsAll[[i]]; if (is.null(E)) return(NULL)
      if (isTRUE(ef$canVanish)) { .fa <- abs(E$f); if (any(!is.finite(.fa)) || min(.fa) < 1e-6 * max(.fa)) return(NULL) }
      nobsAll[i] <- length(E$f)
    }
    totObs <- sum(nobsAll); off <- c(0L, cumsum(nobsAll))
    aB <- matrix(0, totObs, ndir); aReB <- matrix(0, totObs, ndir); aRcB <- matrix(0, totObs, ndir)
    AB <- array(0, c(totObs, ndir, ndir)); R0B <- numeric(totObs); R0sigB <- matrix(0, totObs, nsg)
    dvSensB <- if (length(lamDir)) matrix(0, totObs, ndir) else matrix(0, totObs, 0L)
    jacSum <- setNames(numeric(length(lamNames)), lamNames)
    fB <- numeric(totObs); yB <- numeric(totObs); ehatB <- matrix(0, nsub, neta)
    for (i in seq_len(nsub)) {
      E <- .EsAll[[i]]; E0 <- E0List[[i]]; rows <- (off[i] + 1L):off[i + 1L]
      obs <- .byId[[as.character(.idCode[i])]]; obs <- obs[obs$EVID == 0, , drop = FALSE]
      aB[rows, ] <- E$a; AB[rows, , ] <- E$A; fB[rows] <- E$f; yB[rows] <- .foceiAnalyticTbsY(obs$DV, E$trans)
      if (.fpG) { R0B[rows] <- E$R; aReB[rows, ] <- E$aR; aRcB[rows, ] <- E$aR
        if (nsg > 0L) R0sigB[rows, ] <- E$Rsig }
      else { R0B[rows] <- E0$R; aRcB[rows, ] <- E0$aR                       # aReB stays 0 (frozen)
        if (nsg > 0L) R0sigB[rows, ] <- E0$Rsig }
      if (length(lamDir)) {                            # DV-transform chain (estimated lambda)
        dvSensB[rows, lamDir] <- .foceiAnalyticDvSensLambda(obs$DV, E$trans)
        jacSum <- jacSum + sum(.foceiAnalyticJacLambda(obs$DV, E$trans))
      }
      ehatB[i, ] <- etaSolve[i, ]
    }
    dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
    if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
    ncores <- tryCatch(as.integer(rxode2::getRxThreads()), error = function(e) 1L)
    if (length(ncores) != 1L || is.na(ncores) || ncores < 1L) ncores <- 1L
    .res <- tryCatch(foceiGradAllFoceFR_(aB, AB, aReB, aRcB, R0sigB, dvSensB, fB, yB, R0B, ehatB, as.integer(off),
                                         Oi, dOiCube, if (nom > 0L) as.numeric(tr28) else numeric(0),
                                         neta, nth, nsg, nom, as.integer(dirTh), as.integer(seq_len(nsg)),
                                         as.integer(.fpG), ncores), error = function(e) NULL)
    if (is.null(.res) || !all(is.finite(.res$g)) || !all(is.finite(.res$etaP))) return(NULL)
    g <- .res$g
    for (i in seq_len(nsub)) etaPList[[i]] <- .res$etaP[, , i]
  } else {
    # FOCEI: assemble ALL subjects in ONE OpenMP C++ call (foceiGradAllFR_), removing the
    # per-subject R<->C++ round-trip.  Sensitivities are concatenated over observations.
    nsub <- length(ids); nobsAll <- integer(nsub)
    for (i in seq_len(nsub)) {
      E <- .EsAll[[i]]; if (is.null(E)) return(NULL)
      if (isTRUE(ef$canVanish)) { .fa <- abs(E$f); if (any(!is.finite(.fa)) || min(.fa) < 1e-6 * max(.fa)) return(NULL) }
      nobsAll[i] <- length(E$f)
    }
    totObs <- sum(nobsAll); off <- c(0L, cumsum(nobsAll))       # 0-based per-subject row offsets
    aB <- matrix(0, totObs, ndir); aRB <- matrix(0, totObs, ndir)
    AB <- array(0, c(totObs, ndir, ndir)); ARB <- array(0, c(totObs, ndir, ndir))
    fB <- numeric(totObs); yB <- numeric(totObs); RB <- numeric(totObs)
    RsigB <- matrix(0, totObs, nsg); RsigDirB <- array(0, c(totObs, ndir, nsg))
    dvSensB <- if (length(lamDir)) matrix(0, totObs, ndir) else matrix(0, totObs, 0L)
    jacSum <- setNames(numeric(length(lamNames)), lamNames)
    ehatB <- matrix(0, nsub, neta)
    for (i in seq_len(nsub)) {
      E <- .EsAll[[i]]; rows <- (off[i] + 1L):off[i + 1L]
      obs <- .byId[[as.character(.idCode[i])]]; obs <- obs[obs$EVID == 0, , drop = FALSE]
      aB[rows, ] <- E$a; aRB[rows, ] <- E$aR; AB[rows, , ] <- E$A; ARB[rows, , ] <- E$AR
      fB[rows] <- E$f; yB[rows] <- .foceiAnalyticTbsY(obs$DV, E$trans); RB[rows] <- E$R
      if (nsg > 0L) { RsigB[rows, ] <- E$Rsig; RsigDirB[rows, , ] <- E$RsigDir }
      if (length(lamDir)) {                            # DV-transform chain (estimated lambda)
        dvSensB[rows, lamDir] <- .foceiAnalyticDvSensLambda(obs$DV, E$trans)
        jacSum <- jacSum + sum(.foceiAnalyticJacLambda(obs$DV, E$trans))
      }
      ehatB[i, ] <- etaSolve[i, ]
    }
    dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
    if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
    ncores <- tryCatch(as.integer(rxode2::getRxThreads()), error = function(e) 1L)
    if (length(ncores) != 1L || is.na(ncores) || ncores < 1L) ncores <- 1L
    .res <- tryCatch(foceiGradAllFR_(aB, AB, aRB, ARB, RsigB, RsigDirB, dvSensB, fB, yB, RB, ehatB, as.integer(off),
                                     Oi, dOiCube, if (nom > 0L) as.numeric(tr28) else numeric(0),
                                     neta, nth, nsg, nom, as.integer(dirTh), as.integer(seq_len(nsg)), ncores),
                     error = function(e) NULL)
    if (is.null(.res) || !all(is.finite(.res$g)) || !all(is.finite(.res$etaP))) return(NULL)
    g <- .res$g
    for (i in seq_len(nsub)) etaPList[[i]] <- .res$etaP[, , i]
  }
  names(g) <- c(thStruct, sgNames, omNames)
  if (length(lamNames)) g[lamNames] <- g[lamNames] - 2 * jacSum   # transform Jacobian -2 log|dy'/dDV|
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
  # multiple estimated lambdas (per-endpoint) need an endpoint->lambda DV mapping not yet
  # wired; keep those on FD.  A single estimated lambda is the ported case.
  if (length(.dir$lamNames) > 1L) return(NULL)
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
