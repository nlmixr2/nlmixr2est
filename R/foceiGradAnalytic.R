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

.foceiAnalyticSubjectGradAgqFR <- function(E, Eks, ehat, Om, neta, nth, nsg, dirTh, sigCol,
                                           dOiEst, tr28, ndir, qx, qw, Oi = solve(Om)) {
  tr <- function(M) sum(diag(M))
  .phiU <- function(A) { U <- A; U[lower.tri(U)] <- 0; diag(U) <- diag(U) / 2; U }
  f <- E$f; y <- E$y; R <- E$R; a <- E$a; A <- E$A; aR <- E$aR; AR <- E$AR
  Rsig <- E$Rsig; RsigDir <- E$RsigDir
  res <- y - f
  rf <- -res / R; rR <- 0.5 * (1 / R - res^2 / R^2)
  rff <- 1 / R; rfR <- res / R^2; rRR <- 0.5 * (-1 / R^2 + 2 * res^2 / R^3)
  eff <- 1 / R; eRR <- 0.5 / R^2
  nom <- length(dOiEst); np <- nth + nsg + nom
  ei <- seq_len(neta); di <- seq_len(ndir)
  # exact inner Hessian H (-> etaP) and N: identical to the FOCEI (f,R) kernel
  H <- Oi; N <- matrix(0, neta, ndir)
  for (l in ei) {
    for (m in ei) H[l, m] <- H[l, m] + sum(rff * a[, l] * a[, m] + rfR * (a[, l] * aR[, m] + aR[, l] * a[, m]) +
                                             rRR * aR[, l] * aR[, m] + rf * A[, l, m] + rR * AR[, l, m])
    for (d in di) N[l, d] <- sum(rff * a[, l] * a[, d] + rfR * (a[, l] * aR[, d] + aR[, l] * a[, d]) +
                                   rRR * aR[, l] * aR[, d] + rf * A[, l, d] + rR * AR[, l, d])
  }
  HiM <- solve(H)
  Nsg <- matrix(0, neta, nsg)
  if (nsg > 0L) for (l in ei) for (k in seq_len(nsg)) { .c <- sigCol[k]
    Nsg[l, k] <- sum(rfR * a[, l] * Rsig[, .c] + rRR * aR[, l] * Rsig[, .c] + rR * RsigDir[, l, .c]) }
  # Laplace determinant Hessian Ht == C++ calcEtaHessian H (censOption="gauss":
  # cHff=1/R, cHfr=0, cHrr=0.5/R^2 -- exactly eff/eRR); Ht places the AGQ nodes
  Ht <- Oi; for (l in ei) for (m in ei) Ht[l, m] <- Ht[l, m] + sum(eff * a[, l] * a[, m] + eRR * aR[, l] * aR[, m])
  Hti <- solve(Ht)
  dHtD <- lapply(di, function(s) { D <- matrix(0, neta, neta); for (l in ei) for (m in ei)
    D[l, m] <- sum(-aR[, s] / R^2 * a[, l] * a[, m] + eff * (A[, l, s] * a[, m] + a[, l] * A[, m, s]) +
                     -aR[, s] / R^3 * aR[, l] * aR[, m] + eRR * (AR[, l, s] * aR[, m] + aR[, l] * AR[, m, s])); D })
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
  etaP <- vapply(seq_len(np), function(p) as.numeric(-HiM %*% Mcol(p)), numeric(neta))
  if (neta == 1L) etaP <- matrix(etaP, nrow = 1L)
  # Cholesky factor of Ht and its differential (the only new algebra; dHtStar is
  # the same total derivative the FOCEI trace term uses).  A non-PD/near-singular Ht means
  # the objective placed nodes via the nmNearPD/cholSE branch, not chol(Ht) -- differentiating
  # chol() would be the wrong function, so bail (same PD-margin guard as the batch path).
  H0 <- tryCatch(chol(Ht), error = function(e) NULL)
  if (is.null(H0) || !is.finite(rcond(Ht)) || rcond(Ht) < 1e-10) return(NULL)  # Ht = H0'H0, H0 upper triangular
  Ginv <- backsolve(H0, diag(neta))       # H0^-1
  dHtStarL <- vector("list", np); dGinvL <- vector("list", np)
  for (p in seq_len(np)) {
    .dS <- dHt_p(p); for (l in ei) .dS <- .dS + etaP[l, p] * dHtD[[l]]
    dHtStarL[[p]] <- .dS
    dGinvL[[p]] <- -Ginv %*% .phiU(t(Ginv) %*% .dS %*% Ginv)
  }
  # log joint density up to an eta-independent constant (cancels in pi_k)
  .lf <- function(Ee, eta) -0.5 * sum(log(Ee$R) + (Ee$y - Ee$f)^2 / Ee$R) -
    0.5 * as.numeric(t(eta) %*% Oi %*% eta)
  l0 <- .lf(E, ehat)
  nn <- nrow(qx)
  lognode <- numeric(nn); perNode <- matrix(0, nn, np)
  for (k in seq_len(nn)) {
    x <- qx[k, ]; w <- qw[k, ]; Ek <- Eks[[k]]
    if (is.null(Ek)) return(NULL)
    # sqrt(2): change-of-variable z=sqrt(2)x for the e^{-x^2} Gauss-Hermite kernel, matching
    # inner.cpp's node placement (nodes at sqrt(2)*Ginv*x, untilt exp(x'x)).
    etaCur <- ehat + sqrt(2) * as.numeric(Ginv %*% x)
    lognode[k] <- sum(log(w)) + sum(x^2) + (.lf(Ek, etaCur) - l0)
    .resk <- Ek$y - Ek$f; .Rk <- Ek$R
    rfk <- -.resk / .Rk; rRk <- 0.5 * (1 / .Rk - .resk^2 / .Rk^2)
    # Phi_eta at the NODE: nonzero (the envelope theorem does not apply off the mode)
    gPhik <- as.numeric(Oi %*% etaCur)
    for (l in ei) gPhik[l] <- gPhik[l] + sum(rfk * Ek$a[, l] + rRk * Ek$aR[, l])
    for (p in seq_len(np)) {
      .t <- typ(p)
      .dphi <- if (.t == "th") sum(rfk * Ek$a[, dirTh[p]] + rRk * Ek$aR[, dirTh[p]])
      else if (.t == "sg") sum(rRk * Ek$Rsig[, sigCol[p - nth]])
      # the -tr28 constant is node-independent; sum_k pi_k = 1 passes it through
      else 0.5 * as.numeric(t(etaCur) %*% dOiEst[[omc(p)]] %*% etaCur) - tr28[omc(p)]
      perNode[k, p] <- .dphi + sum(gPhik * (etaP[, p] + sqrt(2) * as.numeric(dGinvL[[p]] %*% x)))
    }
  }
  if (!all(is.finite(lognode))) return(NULL)
  pk <- exp(lognode - max(lognode)); pk <- pk / sum(pk)     # quadrature posterior weights
  g <- vapply(seq_len(np), function(p) 2 * sum(pk * perNode[, p]) + tr(Hti %*% dHtStarL[[p]]),
              numeric(1))
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
                                           ndir, Oi = solve(Om), dvSens = matrix(0, length(E$f), 0L),
                                           censv = integer(0), limv = numeric(0), censOpt = 1L) {
  nom <- length(dOiEst); nobs <- length(E$f)
  dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
  Rsig <- if (is.null(E$Rsig)) matrix(0, nobs, 0L) else E$Rsig
  RsigDir <- if (is.null(E$RsigDir)) array(0, c(nobs, ndir, 0L)) else E$RsigDir
  foceiSubjectGradFR_(E$a, E$A, E$aR, E$AR, Rsig, RsigDir, dvSens,
                      as.integer(censv), as.numeric(limv), as.integer(censOpt), E$f, E$y, E$R,
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
                                               dvSens = matrix(0, length(E$f), 0L),
                                               censv = integer(0), limv = numeric(0)) {
  nobs <- length(E$f); nom <- length(dOiEst)
  .fp <- identical(as.integer(foceType), 1L) || is.null(E0)
  if (.fp) { R0 <- E$R; aRe <- E$aR; aRc <- E$aR; R0sig <- E$Rsig }
  else { R0 <- E0$R; aRe <- matrix(0, nobs, ndir); aRc <- E0$aR; R0sig <- E0$Rsig }
  if (is.null(R0sig)) R0sig <- matrix(0, nobs, nsg)
  dOiCube <- array(0, c(neta, neta, max(nom, 1L)))
  if (nom > 0L) for (k in seq_len(nom)) dOiCube[, , k] <- dOiEst[[k]]
  foceiSubjectGradFoceFR_(E$a, E$A, aRe, aRc, R0sig, dvSens, as.integer(censv), as.numeric(limv),
                          E$f, E$y, R0, as.numeric(ehat), Oi,
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

.foceiOuterFdForFlagged <- function(ids1, analyticRef = matrix(numeric(0), 0, 0)) {
  if (length(ids1) == 0L) return(NULL)
  .g <- tryCatch(foceiOuterFdInd_(as.integer(ids1 - 1L), as.matrix(analyticRef)),
                 error = function(e) NULL)
  if (is.null(.g) || !is.matrix(.g) || nrow(.g) != length(ids1)) return(NULL)
  if (!all(is.finite(.g))) return(NULL)   # a subject that could not be re-optimized
  .g
}

#' Per-FIT constants for the all-C++ analytic outer gradient.
#'
#' Everything the gradient needs that does NOT change between outer iterations: the
#' augmented model's lhs column map, the direction bookkeeping and the problem
#' dimensions.  Computed ONCE per fit and cached in C++, so `analyticOuterGrad` can run
#' with no R interaction at all -- the omega derivatives it also needs come from the
#' fit's own `rxInv` handle, which C++ already holds for the inner problem.
#'
#' Deliberately returns plain atomic vectors: C++ copies them into a POD struct and keeps
#' no R objects alive across the fit.
#' @param ui model UI
#' @param e fit environment (for the live omega/rxInv reuse)
#' @return a plain list, or NULL when the pooled gradient is out of scope
#' @noRd
.foceiGradPooledSetup <- function(ui, e = NULL) {
  tryCatch({
    if (.foceiAnalyticIsMixture(ui)) return(NULL)
    .thv <- tryCatch({
      .ini <- ui$iniDf
      .r <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
      .r <- .r[order(.r$ntheta), , drop = FALSE]
      setNames(as.numeric(.r$est), .r$name)
    }, error = function(.) NULL)
    if (is.null(.thv)) return(NULL)
    ## Omega is only used here for the SHAPE of the estimation-scale derivative block
    ## (how many free omega parameters there are); the values themselves are recomputed
    ## in C++ from the fit's rxInv on every call.  So the initial Omega is fine, and it
    ## is what is available before the fit starts.
    .Om <- tryCatch(get("omega", e), error = function(.) NULL)
    if (is.null(.Om)) .Om <- tryCatch(ui$omega, error = function(.) NULL)
    if (is.null(.Om)) return(NULL)
    st <- .foceiAnalyticGradSetup(ui, .thv, .Om, e)
    if (is.null(st)) return(NULL)
    ## Shape.  `.foceiAnalyticGradSetup` reports interaction = 0 for an ll() model as well
    ## as for FOCE, so isLL has to be carried separately and tested first.
    .isLL <- isTRUE(st$ef$isLL)
    .interaction <- as.integer(st$interaction)
    .nAGQ <- as.integer(st$nAGQ)
    am <- .foceiAnalyticAugModelDirs(ui, st$dir$dirs)
    if (is.null(am)) return(NULL)
    .cols <- tryCatch(.vaeOuterCols(am), error = function(.) NULL)
    if (is.null(.cols)) return(NULL)
    ## The (f,R) kernels contract a variance model.  An ll() endpoint has none -- its
    ## rx_pred_ IS the per-observation log density -- so hasR is required for every shape
    ## except that one.
    if (!.isLL && !isTRUE(.cols$hasR)) return(NULL)
    ## AGQ solves its quadrature nodes through a 1st-order sibling model (26 ODE states
    ## down to 8 on a one-compartment model).  Prefer the one built and disk-cached at
    ## model setup; fall back to building it, and then to the order-2 model, which is
    ## what the nodes used before that optimization existed.
    .colsNode <- NULL
    if (.nAGQ > 1L) {
      .amN <- tryCatch({
        .fmN <- ui$foceiModel
        if (inherits(.fmN$outerNode, "rxode2") && !is.null(.fmN$outerNodeMeta)) {
          c(list(augMod = .fmN$outerNode), .fmN$outerNodeMeta)
        } else {
          .foceiAnalyticAugModelDirs(ui, st$dir$dirs, order = 1L)
        }
      }, error = function(.) NULL)
      if (!is.null(.amN) && identical(.amN$ndir, am$ndir)) {
        .colsNode <- tryCatch(.vaeOuterCols(.amN), error = function(.) NULL)
      }
      if (is.null(.colsNode)) .colsNode <- .cols
    }
    ## Map each OUTER-OPTIMIZER parameter to its slot in the kernel's output vector.
    ##
    ## The kernel emits nth theta directions, then nsg sigma, then nom omega.  That is NOT
    ## the optimizer's parameter vector, for two independent reasons:
    ##   * an estimated boxCox/yeoJohnson lambda appears in BOTH dir$thStruct (it is a
    ##     direction) and ef$sgName (it is a sigma theta of the augmented model), so the
    ##     kernel emits it twice;
    ##   * the mu-referenced (lin/irls) families profile some structural thetas out of the
    ##     outer problem entirely, so the optimizer has FEWER parameters than the kernel.
    ## The deleted R route hid both: it named the vector c(thStruct, sgName, omNames) and
    ## the caller subset it by name, which takes the first match and drops the rest.  C++
    ## is positional and cannot, so carry the map explicitly.  Getting this wrong is
    ## silent -- the arity guard in analyticOuterGradDirect() just declines to finite
    ## differences -- which is how both cases went unnoticed.
    .kernelNames <- c(st$dir$thStruct, st$ef$sgName, st$omNames)
    .thAll <- ui$iniDf$name[!is.na(ui$iniDf$ntheta)]
    .thRows <- ui$iniDf[!is.na(ui$iniDf$ntheta), , drop = FALSE]
    .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
    .thEst <- .thRows$name[!.thRows$fix]
    .parNames <- c(setdiff(.thEst, .foceiMuSkipThetaNames(ui, .thAll)), st$omNames)
    .gMap <- match(.parNames, .kernelNames)
    if (anyNA(.gMap)) return(NULL)
    .gMap <- as.integer(.gMap - 1L)          # 0-based for C++
    ## ntheta position of each structural theta -- the ll() perturbation of a non-mu
    ## theta moves th[thPos[p]], which is not the direction index.
    .thPos <- tryCatch(as.integer(ui$iniDf$ntheta[match(st$dir$thStruct, ui$iniDf$name)]),
                       error = function(.) integer(0))
    list(cols = .cols,
         colsNode = .colsNode,
         neta = as.integer(st$neta), nth = as.integer(st$dir$nth),
         nsg = as.integer(length(st$ef$sgName)), nom = as.integer(length(st$dOiEst)),
         dirTh = as.integer(st$dir$dirTh),
         sigCol = as.integer(seq_len(length(st$ef$sgName))),
         lamDir = as.integer(st$dir$lamDir),
         nLam = as.integer(length(st$dir$lamNames)),
         censOpt = as.integer(rxode2::rxGetControl(ui, "censOption", 0L)),
         isLL = .isLL,
         interaction = .interaction,
         foceType = as.integer(st$foceType),
         nAGQ = .nAGQ,
         ebeTol = {
           ## FOCE frozen-R0 Newton score tolerance -- a fixed 1e-9, NOT derived from
           ## sigdig.  foceiControl(foceEbeTol=) overrides it.
           .t <- suppressWarnings(as.numeric(rxode2::rxGetControl(ui, "foceEbeTol", NA_real_)))
           if (!is.finite(.t) || .t <= 0) .t <- 1e-9
           .t
         },
         ebeSkipTol = 1e-3,   ## looser first-iteration "already stationary?" test
         dependsF0 = isTRUE(st$ef$dependsF0),
         canVanish = isTRUE(st$ef$canVanish),
         thPos = .thPos,
         gMap = .gMap)
  }, error = function(e) NULL)
}

#' The FIT's own ODE tolerances, for solving the analytic gradient.
#'
#' There is deliberately no separate "analytic gradient" tolerance: the gradient has to
#' be the gradient of the objective the optimizer is minimizing, so the augmented model
#' is solved at the tolerance that objective is solved at.  Tightening it here makes the
#' gradient describe a different objective -- and the finite-difference references these
#' gradients are checked against are themselves objective differences at the fit's
#' tolerance.  [.foceiAnalyticSolveTol()] remains correct for the COVARIANCE path, where
#' `covSolveTol` is a documented user control.
#' @return `c(atol, rtol)`, or NULL when the fit's values cannot be read
#' @noRd
.foceiGradSolveTol <- function(ui) {
  .rc <- tryCatch(ui$control$rxControl, error = function(e) NULL)
  if (is.null(.rc)) return(NULL)
  .a <- suppressWarnings(as.numeric(.rc$atol)[1])
  .r <- suppressWarnings(as.numeric(.rc$rtol)[1])
  if (!is.finite(.a) || !is.finite(.r) || .a <= 0 || .r <= 0) return(NULL)
  c(.a, .r)
}

.foceiGradSolveTolOr <- function(ui) {
  .t <- .foceiGradSolveTol(ui)
  if (is.null(.t)) .foceiAnalyticSolveTol(ui) else .t
}

#' Subjects whose augmented solve failed, for the Phase 8D2 per-individual FD.
#'
#' Recorded here rather than acted on in the solve loop: that loop runs inside
#' OdeSwapEsBatch(odeSlotOuter), and the finite difference needs the INNER problem's
#' event sensitivities, which can only be installed at a batch boundary.
#' @noRd
.foceiOuterFlagged <- new.env(parent = emptyenv())
.foceiOuterFlagged$ids <- integer(0)

.foceiAnalyticSolveAll <- function(am, thv, ebes, ids, data, obsTimes, tol = 1e-10) {
  ## Solve the augmented model IN THE SHARED FOCEi pool (which it sized) and take
  ## the per-subject E structures straight from C++, instead of routing through
  ## rxode2::rxSolve, which frees and rebuilds the global solve on every call.
  ## Used by BOTH est="vae"'s M-step and focei's own fast gradient.
  ##
  ## No session flag guards this any more.  vaeOuterSolve_ refuses unless the
  ## augmented model is registered AND the pool is at least its size
  ## (odeSwapCanPool -> odeDenyPoolNotSized otherwise), which is the structural
  ## form of what `.vaeGradEnv$active` was patching: a focei fast fit after a vae
  ## grad fit, running against a pool sized for its own inner model.  A NULL
  ## falls through to the rxSolve path below.
  ##
  ## DDE is IN scope for the pooled route.  The old exclusion assumed a delay model
  ## pins method="dop853"/dense per solve, which a shared pool cannot do -- but focei
  ## already forces that configuration at the FIT level (R/focei.R, the hasDelay block:
  ## method 0, stiff2 13, dense TRUE), so a DDE fit's pool is built that way to begin
  ## with and there is nothing to change per solve.
  ## (This used to be gated by .odeSwapNoPool, a verification-only opt-out that let a
  ## test evaluate the same fit through rxSolve instead of the pool.  Its only setter was
  ## the R gradient route, which is gone, so the gate could never fire.)
  {
    .cols <- tryCatch(.vaeOuterCols(am), error = function(e) NULL)
    if (!is.null(.cols)) {
      .nc <- tryCatch({ .c <- am$cores
        if (is.null(.c) || is.na(.c) || .c < 1L) 1L else as.integer(.c) },
        error = function(e) 1L)
      ## No tolerance argument: the pooled solve runs at the FIT's tolerance (C++
      ## OdeFitTolGuard resets to it), because the gradient must be the gradient of the
      ## objective the optimizer is minimizing.  `tol` below still tunes the rxSolve
      ## fallback and the covariance path.
      .Ec <- tryCatch(vaeOuterSolve_(as.numeric(thv), as.matrix(ebes), .cols, .nc),
                      error = function(e) NULL)
      ## vaeOuterSolve_ now flags failed subjects per individual (attr "ok") rather
      ## than discarding the whole gradient.  Until the Phase 8D2 finite-difference
      ## phase consumes those flags, a flagged subject still falls through to the
      ## rxSolve route, i.e. behaviour is unchanged -- but the flags are here.
      if (!is.null(.Ec) && length(.Ec) > 0L) {
        .ok <- attr(.Ec, "ok")
        if (is.null(.ok) || all(.ok == 1L)) {
          .foceiOuterFlagged$ids <- integer(0)
          return(.Ec)
        }
        ## A flagged subject has no E.  Give it a zero-filled one of the right shape so
        ## the assembly below keeps working; its gradient column is replaced wholesale in
        ## foceiGradAllFR_ by the per-individual finite difference, so these zeros never
        ## reach the result.
        .good <- which(.ok == 1L)
        if (length(.good) > 0L) {
          .tmpl <- .Ec[[.good[1L]]]
          for (.i in which(.ok == 0L)) {
            .nobsI <- length(obsTimes[[.i]])
            .z <- lapply(.tmpl, function(.x) {
              if (is.null(dim(.x))) numeric(.nobsI)
              else array(0, c(.nobsI, dim(.x)[-1]))
            })
            .z <- .z[names(.tmpl)]
            if (!is.null(.tmpl$trans)) .z$trans <- .tmpl$trans
            .Ec[[.i]] <- .z
          }
          .foceiOuterFlagged$ids <- which(.ok == 0L)
          attr(.Ec, "ok") <- .ok
          return(.Ec)
        }
      }
    }
  }
  dirs <- am$dirs; nd <- length(dirs); neta <- ncol(ebes)
  etav <- paste0("ETA_", seq_len(neta), "_")
  pars <- data.frame(ID = ids)
  for (k in seq_len(neta)) pars[[etav[k]]] <- ebes[, k]
  for (.nm in names(thv)) pars[[.nm]] <- thv[[.nm]]
  .ev <- .foceiAnalyticEvents(am, data)                    # reuse the pre-translated event table
  .nc <- if (is.null(am$cores)) 0L else am$cores           # fit's rxControl thread count (parallel)
  # DDE: force pure dop853 (dense, no Jacobian) -- its 8th-order dense history reproduces the
  # delayed sensitivity solve exactly and needs no Jacobian, sidestepping the composite/ros4
  # on-the-fly Jacobian generation for this THETA/ETA-named augmented model.
  .ddeArgs <- if (isTRUE(rxode2::rxModelVars(am$augMod)$flags[["hasDelay"]] == 1L))
    list(method = "dop853", stiff2 = 0L, dense = TRUE) else list()
  .sol <- tryCatch(withCallingHandlers(
    as.data.frame(do.call(rxode2::rxSolve, c(list(am$augMod, params = pars, events = .ev, cores = .nc,
                                  returnType = "data.frame", atol = tol[1], rtol = tol[length(tol)]), .ddeArgs))),
    warning = function(w) invokeRestart("muffleWarning")), error = function(e) NULL)
  if (is.null(.sol) || !all(c("rx_predf_", paste0("rx_f1_", if (is.null(am$fDirs)) dirs else am$fDirs)) %in% names(.sol))) return(NULL)
  # Extract every sensitivity column from the WHOLE solve as a matrix ONCE (the per-subject
  # data.frame [[ + paste0 dominated the gradient -- the ODE solve itself is ~4%); slice by row
  # per subject.  Column names + index maps are precomputed on `am$cols` at build time.
  .cm <- if (is.null(am$cols)) .foceiAnalyticCols(dirs, if (is.null(am$fDirs)) dirs else am$fDirs, am$P2, if (is.null(am$P2r)) am$P2 else am$P2r, am$sigTh) else am$cols
  np2 <- nrow(am$P2); np2r <- if (is.null(am$P2r)) np2 else nrow(am$P2r)
  .hasR <- isTRUE(am$hasRvar); .hasT <- isTRUE(am$hasTrans)
  # Every column subset below must be present, or `.sol[, cols]` throws "undefined columns
  # selected" (e.g. an rxode2 version/feature that did not emit an rvar/rsig/transform column).
  # Verify up front and cleanly return NULL so the caller falls back to finite differences.
  .need <- c("rx_predf_", "time", .cm$f1, .cm$f2)
  if (.hasR) { .need <- c(.need, "rx_rvarf_", .cm$rvar1, .cm$rvar2)
    if (length(am$sigTh) > 0L) .need <- c(.need, .cm$rsig, unlist(.cm$rsig1), .cm$rsig2) }
  if (.hasT) .need <- c(.need, "rx_tyj_", "rx_tlambda_", "rx_tlow_", "rx_thi_")
  if (!all(.need %in% names(.sol))) return(NULL)
  .M1 <- as.matrix(.sol[, .cm$f1, drop = FALSE]); .M2 <- as.matrix(.sol[, .cm$f2, drop = FALSE])
  .fp <- .sol$rx_predf_; .tm <- .sol$time
  if (.hasR) {
    .MR1 <- as.matrix(.sol[, .cm$rvar1, drop = FALSE]); .MR2 <- as.matrix(.sol[, .cm$rvar2, drop = FALSE])
    .Rf <- .sol$rx_rvarf_; .nsig <- length(am$sigTh)
    if (.nsig > 0L) { .MRs <- as.matrix(.sol[, .cm$rsig, drop = FALSE])
      .MRs1 <- lapply(.cm$rsig1, function(cc) as.matrix(.sol[, cc, drop = FALSE]))
      .MRs2 <- as.matrix(.sol[, .cm$rsig2, drop = FALSE]) }
  }
  if (.hasT) .tr <- list(yj = .sol$rx_tyj_, lambda = .sol$rx_tlambda_, low = .sol$rx_tlow_, hi = .sol$rx_thi_)
  .idcol <- if ("id" %in% names(.sol)) .sol$id else .sol$ID
  .byIdSol <- split(seq_len(nrow(.sol)), as.character(.idcol))
  Es <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    .ri <- .byIdSol[[as.character(ids[i])]]
    if (is.null(.ri)) .ri <- .byIdSol[[as.character(i)]]
    if (is.null(.ri)) return(NULL)
    .keep <- .ri[.tm[.ri] %in% obsTimes[[i]]]
    no <- length(.keep); if (no != length(obsTimes[[i]])) return(NULL)
    a <- matrix(0, no, nd); a[, .cm$fDirIdx] <- .M1[.keep, , drop = FALSE]; A <- array(0, c(no, nd, nd))
    for (r in seq_len(np2)) { .v <- .M2[.keep, r]; A[, .cm$iiF[r], .cm$jjF[r]] <- .v; A[, .cm$jjF[r], .cm$iiF[r]] <- .v }
    .E <- list(f = .fp[.keep], a = a, A = A)
    if (.hasR) {
      aR <- .MR1[.keep, , drop = FALSE]; AR <- array(0, c(no, nd, nd))
      for (r in seq_len(np2r)) { .v <- .MR2[.keep, r]; AR[, .cm$ii[r], .cm$jj[r]] <- .v; AR[, .cm$jj[r], .cm$ii[r]] <- .v }
      .E$R <- .Rf[.keep]; .E$aR <- aR; .E$AR <- AR
      if (.nsig > 0L) {
        .E$Rsig <- .MRs[.keep, , drop = FALSE]
        .E$RsigDir <- array(vapply(.MRs1, function(M) M[.keep, , drop = FALSE], matrix(0, no, nd)), c(no, nd, .nsig))
        .Rs2 <- array(0, c(no, .nsig, .nsig))
        for (r in seq_len(nrow(.cm$sigP2))) { .v <- .MRs2[.keep, r]
          .Rs2[, .cm$sigP2$a[r], .cm$sigP2$b[r]] <- .v; .Rs2[, .cm$sigP2$b[r], .cm$sigP2$a[r]] <- .v }
        .E$Rsig2 <- .Rs2
      }
    }
    if (.hasT) .E$trans <- lapply(.tr, `[`, .keep)       # both-sides transform: DV -> tbs(DV) scale
    Es[[i]] <- .E
  }
  Es
}

.foceiAnalyticIsMixture <- function(ui) {
  isTRUE(tryCatch(length(ui$thetaMixIndex) > 0L, error = function(e) FALSE))
}

.foceiAnalyticGradSetup <- function(ui, thVals, Om, e = NULL,
                                    caller = .analyticGradCaller(ui)) {
  if (is.na(caller)) return(NULL)
  if (!.hasRxSens()) return(NULL)
  if (.foceiAnalyticIsMixture(ui)) return(NULL)     # mixtures: weighted sum, no treatment yet
  if (isTRUE(any(ui$predDf$linCmt))) return(NULL)   # linCmt(): no symbolic state sensitivities
  if (!.analyticGradAllowsBoundedTr(ui, caller)) return(NULL)
  # tad/podo/tafd/tlast/tfirst/dosenum are functions of time and the dose record
  # only (no eta/theta dependence), so rxode2 treats them as zero-derivative
  # constants in the sensitivity expansion (.rxToSEDualVarFunction) -- they no
  # longer need to force the finite-difference fallback.
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(NULL)
  # ll()/generalized likelihood (needOptimHess -> interaction=0, EXACT inner
  # Hessian): rx_pred_ is the log-density, so skip the Gaussian ErrFull and set up
  # the direct-log-density core (gradPooledCoreLL in C++).
  if (.foceiLLGradInScope(ui, caller)) {
    .map <- .foceiEtaThetaMap(ui); neta <- length(.map$etaNames)
    .dir <- .foceiOuterDirsLL(ui); if (is.null(.dir)) return(NULL)
    .oe <- .foceiEstOmegaDeriv(ui, Om, e); if (is.null(.oe)) return(NULL)
    return(list(ef = list(isLL = TRUE), dir = .dir, dOiEst = .oe$dOi, tr28 = .oe$tr28,
                omNames = .oe$names, neta = neta, etaNames = .map$etaNames,
                interaction = 0L, foceType = 0L, nAGQ = 1L))
  }
  interaction <- as.integer(rxode2::rxGetControl(ui, "interaction", 1L))            # 1 FOCEI / 0 FOCE
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(ui, "foceType", 0L)) else 0L
  ## FOCE (interaction = 0) was declined here up front: its frozen-R0 EBE Newton could
  ## not reach the 1e-9 score target at the default solve, |S| flooring near 5e-3 at
  ## rtol = 1e-3 (nlmixr2/nlmixr2est#836).  That was measured BEFORE the shared ODE solve
  ## pool was fixed (#839), where a peer solve run under another slot's event-sensitivity
  ## shape corrupted the scratch the score is assembled from.  The gate is lifted so FOCE
  ## goes through gradPooledCore's isFoce/foceEbeNewton path like any other shape; a
  ## Newton that still cannot converge declines per fit at its own site rather than
  ## being refused for the whole method.
  nAGQ <- as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L))
  # agqControl() forces interaction=TRUE, so only the FOCEI (f,R) kernel has a quadrature
  # form -- a FOCE-AGQ combination cannot arise.
  if (nAGQ > 1L && interaction != 1L) return(NULL)
  # the aqLow/aqHi clamp (inner.cpp) kinks the objective; both default to +/-Inf
  if (nAGQ > 1L &&
        (is.finite(as.numeric(rxode2::rxGetControl(ui, "agqLow", -Inf))) ||
           is.finite(as.numeric(rxode2::rxGetControl(ui, "agqHi", Inf))))) return(NULL)
  # The grid is placed by Ht's Cholesky FACTOR (etaCur = etahat + chol(Ht)^-1 x), so we
  # must differentiate the exact factorization the objective used.  cholSEOpt forces the
  # generalized Cholesky, whose factor differs from chol() even for a PD matrix.  (The
  # runtime doChol flips are temporary and fire only when the plain Cholesky fails --
  # where chol(Ht) fails in R too, so those are already covered.)
  if (nAGQ > 1L && isTRUE(as.logical(rxode2::rxGetControl(ui, "cholSEOpt", FALSE)))) return(NULL)
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
  .dir <- .foceiAnalyticDirections(ini, thetaForEta, ef$sgName, neta,
                                   sharedEta = unname(.foceiEtaOccurrence(ui) > 1L))
  if (is.null(.dir)) return(NULL)
  # multiple estimated lambdas (per-endpoint) need an endpoint->lambda DV mapping not yet
  # wired; keep those on FD.  A single estimated lambda is the ported case.
  if (length(.dir$lamNames) > 1L) return(NULL)
  .oe <- .foceiEstOmegaDeriv(ui, Om, e); if (is.null(.oe)) return(NULL)
  list(ef = ef, dir = .dir, dOiEst = .oe$dOi, tr28 = .oe$tr28, omNames = .oe$names,
       neta = neta, etaNames = etaNames, interaction = interaction, foceType = foceType,
       nAGQ = nAGQ)
}

#' Post-fit analytic natural-scale gradient, computed by the fit's OWN C++ path.
#'
#' Re-enters the estimation machinery at the fit's converged estimates and EBEs --
#' `est="none"` with zero inner/outer iterations and the fit's `etaMat`, the same
#' post-fit re-entry `setCov()` uses -- and reads back the gradient
#' `analyticOuterGrad()` stashed.  So this is the SHIPPING gradient, not a parallel
#' implementation of it: that distinction is the whole point, because a test that
#' validates a second implementation proves nothing about the code the fit runs.
#'
#' `NULL` when the fit is out of analytic scope -- the same signal
#' `foceiControl(fast=TRUE)` acts on when it falls back to finite differences -- and
#' also when the fit was not run with `fast=TRUE` at all.
#' @param fit nlmixr2 fit object
#' @return named natural-scale gradient (structural thetas, sigmas, om.chol), or `NULL`
#' @noRd
.foceiGradDirect <- function(fit) {
  tryCatch({
    .env <- if (rxode2::rxIs(fit, "nlmixr2FitData")) fit$env else fit
    .est <- .env$est
    if (is.null(.est) || !nzchar(.est)) return(NULL)
    .control <- .env$foceiControl
    .control$maxInnerIterations <- 0L      # evaluate at the fit's EBEs, do not re-optimize
    .control$maxOuterIterations <- 0L      # no outer step: the gradient is at THIS theta
    .control$calcTables <- FALSE
    .control$covMethod <- 0L               # no covariance step
    .control$skipCov <- fit$skipCov
    # `fast` is deliberately NOT forced on: this reports what the analytic gradient does
    # for THIS fit as configured, so a fast=FALSE fit correctly yields NULL rather than a
    # gradient it never used.
    # Re-run under the fit's OWN est, not est="none": the gradient SHAPE is the
    # estimation method (FOCE freezes the residual variance, AGQ adds quadrature), and
    # est="none" would silently evaluate every fit as plain FOCEI.
    .ui <- fit$ui
    .th <- tryCatch(fit$theta, error = function(.) NULL)
    if (!is.null(.th)) {                   # pin the final thetas on the ui
      .w <- match(names(.th), .ui$iniDf$name); .ok <- !is.na(.w)
      .ui$iniDf$est[.w[.ok]] <- as.numeric(.th)[.ok]
    }
    .eta <- tryCatch(fit$eta, error = function(.) NULL)
    if (!is.null(.eta)) {                  # ...and the final etas
      .control$etaMat <- as.matrix(.eta[, setdiff(names(.eta), "ID"), drop = FALSE])
    }
    # the nested re-fit resets mu-referencing global state (.muRefTrans$cur); restore it
    .savedMuRef <- .muRefTrans$cur
    on.exit(.muRefTrans$cur <- .savedMuRef, add = TRUE)
    .f2 <- suppressMessages(suppressWarnings(
      nlmixr2(.ui, data = getData(fit), est = .est, control = .control)))
    .src <- tryCatch(.f2$env, error = function(.) NULL)
    if (is.null(.src) || !exists(".gradDirectFirst", .src, inherits = FALSE)) return(NULL)
    .g <- as.numeric(get(".gradDirectFirst", .src))
    # Name it the way the gradient assembly orders it: structural thetas, then sigmas,
    # then the estimation-scale omega (Cholesky) elements.
    .ini <- .ui$iniDf
    .thRows <- .ini[!is.na(.ini$ntheta), , drop = FALSE]
    .thRows <- .thRows[order(.thRows$ntheta), , drop = FALSE]
    .thv <- fit$theta[.thRows$name]
    if (anyNA(.thv)) .thv <- .thRows$est
    .st <- .foceiAnalyticGradSetup(.ui, stats::setNames(as.numeric(.thv), .thRows$name),
                                   fit$omega)
    if (is.null(.st)) return(NULL)
    # These names are in KERNEL space (nth + nsg + nom).  That is not the outer
    # optimizer's vector whenever a parameter occupies two kernel slots -- an estimated
    # boxCox/yeoJohnson lambda is both a theta direction and a sigma slot, so the kernel
    # names come out one longer than the gradient (9 vs 8 on a 1-cmt boxCox model) and
    # this used to bail to NULL, reporting "no analytic gradient" for a gradient that
    # had in fact been computed.  gMap is the same kernel -> outer gather the C++ uses
    # (analyticOuterGradDirect), so reuse it rather than re-deriving the correspondence.
    .nmKer <- c(.st$dir$thStruct, .st$ef$sgName, .st$omNames)
    .nm <- .nmKer
    if (length(.nmKer) != length(.g)) {
      .gp <- tryCatch(.foceiGradPooledSetup(.ui), error = function(e) NULL)
      .map <- if (is.null(.gp)) NULL else .gp$gMap
      if (is.null(.map) || length(.map) != length(.g) ||
            any(.map < 0L) || any(.map >= length(.nmKer))) return(NULL)
      .nm <- .nmKer[.map + 1L]
    }
    stats::setNames(.g, .nm)
  }, error = function(e) NULL)
}

.analyticGradCaller <- function(ui) {
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fast", FALSE)))) return("focei")
  if (identical(as.character(rxode2::rxGetControl(ui, "nonMuTheta", "")), "grad")) return("vae")
  NA_character_
}

#' Bounded-transform scope gate.
#'
#' `preProcessBoundedTransform` records the transforms on the ALREADY-REWRITTEN
#' ui, so by the time the gradient sees them the model is on the unconstrained
#' `rxBoundedTr.*` scale.  focei must still bail: it REPORTS a natural-scale
#' gradient to the outer optimizer, which would need a Jacobian correction that is
#' not applied.  The VAE consumes the gradient internally, on the same
#' unconstrained scale it takes its M-step on, so no correction arises.
#' @noRd
.analyticGradAllowsBoundedTr <- function(ui, caller) {
  if (identical(caller, "vae")) return(TRUE)
  is.null(ui$boundedTransforms) || length(ui$boundedTransforms) == 0L
}

#' Is the analytic outer gradient in scope for a VAE fit?
#'
#' Cheap direction-set probe -- no symengine/gcc pass -- covering every static
#' gate: `linCmt()`, `fo`, the distribution/error-model scope, IOV, and a model
#' with no eta.  A later build or solve failure still falls back at runtime.
#'
#' Two admissible shapes, the same pair `.foceiAnalyticGradSetup` dispatches on:
#' a conditionally Gaussian endpoint (the `(f,R)` direction set) or a single
#' non-Gaussian `ll()`/generalized endpoint (the direct-log-density set).  The
#' VAE consumes either through the SAME C++ gradient core, which
#' routes on `ef$isLL`.
#' @noRd
.vaeGradInScope <- function(ui) {
  if (!is.null(tryCatch(.foceiOuterDirs(ui, "vae"), error = function(e) NULL))) return(TRUE)
  isTRUE(.foceiLLGradInScope(ui, "vae")) &&
    !is.null(tryCatch(.foceiOuterDirsLL(ui), error = function(e) NULL))
}

#' Direction set for the augmented outer-gradient model, computed from the UI
#' alone (does not depend on theta/eta values): one direction per eta plus one per
#' non-mu-referenced structural theta.  `NULL` if out of analytic scope.
#' @noRd
.foceiOuterDirs <- function(ui, caller = .analyticGradCaller(ui)) {
  if (!.hasRxSens()) return(NULL)
  if (.foceiAnalyticIsMixture(ui)) return(NULL)     # mixtures: weighted sum, no treatment yet
  if (isTRUE(any(ui$predDf$linCmt))) return(NULL)   # linCmt(): no symbolic state sensitivities
  if (!.analyticGradAllowsBoundedTr(ui, caller)) return(NULL)
  if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(NULL)
  ef <- .foceiAnalyticErrFull(ui); if (is.null(ef)) return(NULL)
  .map <- .foceiEtaThetaMap(ui); neta <- length(.map$etaNames)
  if (neta == 0L) return(NULL)
  if (length(.uiIovEnv$iovVars) > 0L) return(NULL)
  .foceiAnalyticDirections(ui$iniDf, .map$thetaForEta, ef$sgName, neta,
                           sharedEta = unname(.foceiEtaOccurrence(ui) > 1L))
}

#' Is a fit in scope for the ll()/generalized-likelihood analytic outer gradient?
#' The ll() objective uses the EXACT inner Hessian (needOptimHess), so `rx_pred_`
#' is the per-observation log-density and the gradient is assembled by
#' the log-density core (differentiating it directly) rather than
#' the Gaussian (f,R) path.  Phase 1 scope: a single non-Gaussian endpoint, no
#' linCmt/bounded transform/IOV/FO, at least one eta.  (Multi-endpoint, censoring,
#' and nAGQ are handled by falling back to the finite-difference gradient.)
#'
#' `caller` only reaches the bounded-transform gate, which focei must fail and
#' the VAE need not -- see `.analyticGradAllowsBoundedTr`.  Defaulted, so the
#' focei callers keep their exact behavior.
#' @noRd
.foceiLLGradInScope <- function(ui, caller = .analyticGradCaller(ui)) {
  tryCatch({
    if (!.hasRxSens()) return(FALSE)
    .pd <- ui$predDfFocei
    if (is.null(.pd) || nrow(.pd) != 1L) return(FALSE)
    ## Single endpoint.  Multi-endpoint ll() was TRIED and reverted: lifting this gate
    ## yields a gradient that runs (nAnalyticGradDirect > 0, "grad: analytic") but does
    ## NOT verify -- on a 2-endpoint warfarin ll() model, analytic vs central differences
    ## was off by 7.6x on tcl and ~370x on add.pd.  So the mechanism works and the
    ## mathematics does not.  Unlike the GAUSSIAN multi-endpoint case (which the R route
    ## already served, and which the pooled route reproduces bit-identically), there is
    ## no shipping reference here to fall back on.
    ##
    ## Re-measured after the solve-pool fix (#839) and the multi-endpoint pooling / CMT
    ## re-basing: tcl is now correct (7.6x -> 0.36%), but add.pd is still ~373x off and
    ## tka/tv/add.pk are 4.2x/1.9x/2.5x off.  What is right is tcl (the only structural
    ## theta carrying an eta) plus the PD-only algebraic thetas; what is wrong is the
    ## non-eta structural thetas and the per-endpoint RESIDUAL thetas.  That is direction
    ## bookkeeping, not pooling -- see the thPos/gMap note in .foceiGradPooledSetup, which
    ## is single-endpoint-shaped.  nlmixr2/nlmixr2est#838.
    if (all(as.character(.pd$distribution) %in% c("norm", "dnorm"))) return(FALSE)  # Gaussian -> (f,R) path
    # loadPruneSens clears predDfFocei$linCmt for a promoted solved-form linCmt(), so it
    # passes this coarse scope gate.  Its 1st-order eta sensitivity converts (rxode2
    # linCmtB), but the 2nd-order does NOT (rxFromSE cannot emit the nested linCmtB
    # derivative), so .foceiAddHdEta2 fails and the fit falls back to the finite-difference
    # Hessian/gradient at build time (see .foceiMaybeAddHdEta2).  A residual TRUE here marks
    # a case the promotion cannot cover -- out of scope like the Gaussian path.
    if (isTRUE(any(ui$predDfFocei$linCmt))) return(FALSE)
    if (!.analyticGradAllowsBoundedTr(ui, caller)) return(FALSE)
    if (isTRUE(as.logical(rxode2::rxGetControl(ui, "fo", FALSE)))) return(FALSE)
    if (as.integer(rxode2::rxGetControl(ui, "nAGQ", 1L)) > 1L) return(FALSE)
    if (length(.uiIovEnv$iovVars) > 0L) return(FALSE)
    length(.foceiEtaThetaMap(ui)$etaNames) > 0L
  }, error = function(e) FALSE)
}

#' Direction set for the ll() analytic outer gradient: one direction per eta plus
#' one per non-mu-referenced structural theta.  For an ll() endpoint there is no
#' Gaussian residual-sigma set (add.sd et al. appear directly in the log-density),
#' so every non-mu structural theta gets its own direction (`sgName = character(0)`).
#' @noRd
.foceiOuterDirsLL <- function(ui) {
  .map <- .foceiEtaThetaMap(ui); neta <- length(.map$etaNames)
  if (neta == 0L) return(NULL)
  .d <- .foceiAnalyticDirections(ui$iniDf, .map$thetaForEta, character(0), neta,
                                 sharedEta = unname(.foceiEtaOccurrence(ui) > 1L))
  if (is.null(.d) || is.null(.d$dirs)) return(NULL)
  .d
}

# Build the augmented outer-gradient sensitivity model (compiled model + `dirs` +
# `P2`) for a UI.  This is the persistent `..outer` sibling of the inner model:
# it depends only on the model + direction set (NOT theta/eta/omega), so it is
# built once during model setup (via `rxUiGet.foceiModel`/`foceModel`, which
# disk-cache the whole model list) and reused across every outer-gradient call.
# Callable independently as `ui$foceiOuter`.  `NULL` when out of analytic scope
# (the gradient then falls back to finite differences).
#' @export
rxUiGet.foceiOuter <- function(x, ...) {
  .ui <- x[[1]]
  .caller <- .analyticGradCaller(.ui)
  if (is.na(.caller)) return(NULL)
  interaction <- as.integer(rxode2::rxGetControl(.ui, "interaction", 1L))
  foceType <- if (interaction == 0L) as.integer(rxode2::rxGetControl(.ui, "foceType", 0L)) else 0L
  # nAGQ > 1 (adaptive Gaussian quadrature) uses the SAME augmented model at eta-hat: the
  # quadrature nodes are extra eta points on the same sensitivity solve, so the
  # direction set and the symbolic expansion are unchanged.  (The nodes themselves solve
  # a cheaper 1st-order model -- see rxUiGet.foceiOuterNode.)
  .dir <- .foceiOuterDirs(.ui, .caller)
  # ll()/generalized endpoint (needOptimHess, interaction=0): the Gaussian (f,R)
  # direction builder declines (ErrFull is norm-only), but rx_pred_ is the
  # log-density and the same augmented model supplies its 1st/2nd-order eta/theta
  # derivatives -- build over the ll() direction set instead.
  if (is.null(.dir) && .foceiLLGradInScope(.ui, .caller)) .dir <- .foceiOuterDirsLL(.ui)
  if (is.null(.dir)) return(NULL)
  .foceiAnalyticAugModelDirs(.ui, .dir$dirs)
}
attr(rxUiGet.foceiOuter, "rstudio") <- emptyenv()

#' Augmented model for the AGQ quadrature NODES: the same direction set as `foceiOuter`
#' but 1st order only.
#'
#' The nodes read only `f`/`R`/`a`/`aR`/`Rsig` -- they never touch `A`/`AR`/`RsigDir`,
#' because the 2nd-order block is used solely at eta-hat (the exact inner Hessian -> etaP,
#' and dHtD).  Dropping that tier takes the node solve from 26 ODE states to 8 on a
#' one-compartment model, and the nodes are `nAGQ^neta` solves per gradient -- 45-54% of
#' the gradient at neta=3/nAGQ=3 and 77-86% at neta=5.  Measured: 1.07x (neta=3, nAGQ=2)
#' to 1.91x (neta=5, nAGQ=3) on the whole gradient.
#'
#' Only built for nAGQ > 1; every other fast fit gets NULL and pays no extra build.  Like
#' `foceiOuter` this rides in the disk-cached `foceiModel` list, so the extra symengine+gcc
#' pass is paid once per model, not once per session.
#' @noRd
#' @export
rxUiGet.foceiOuterNode <- function(x, ...) {
  .ui <- x[[1]]
  if (!isTRUE(rxode2::rxGetControl(.ui, "fast", FALSE))) return(NULL)
  if (as.integer(rxode2::rxGetControl(.ui, "nAGQ", 1L)) <= 1L) return(NULL)
  .dir <- .foceiOuterDirs(.ui); if (is.null(.dir)) return(NULL)
  .foceiAnalyticAugModelDirs(.ui, .dir$dirs, order = 1L)
}
attr(rxUiGet.foceiOuterNode, "rstudio") <- emptyenv()

#' Estimation-scale (Cholesky) Omega-inverse derivatives for the outer gradient's
#' Omega block, from rxode2's `rxSymInvCholEnvCalculate` (rxSymInv.R d.omegaInv /
#' tr.28).  Returns `list(dOi = list(dOmega^{-1}/d theta_omega_k), tr28 = 0.5 *
#' tr(dOmega^{-1}_k Omega), names = <omega parameter names>)`, or `NULL`.
#' @noRd
.foceiEstOmegaDeriv <- function(ui, Om, e = NULL) {
  tryCatch({
    # Build the rxSymInvChol env from the current Omega with the SAME diagonal
    # transform the optimizer uses, so the Cholesky parameter order/scale matches
    # op_focei's Omega slots (rxUiGet.focei builds env$rxInv the same way).
    .diagXform <- rxode2::rxGetControl(ui, "diagXform", "sqrt")
    # A fresh env pays ~60ms of one-time symbolic setup on its first `$d.omegaInv` (a
    # repeat read is free), which made this ~40% of each analytic gradient.  Reuse the
    # fit's persistent `env$rxInv` (C++ keeps it current via setOmegaTheta) -- but only
    # when it is already at this Omega, so a stale env falls back instead of silently
    # returning derivatives at the wrong one.
    .rxInv <- NULL
    if (!is.null(e)) {
      .cand <- tryCatch(get("rxInv", e), error = function(.) NULL)
      if (!is.null(.cand) && rxode2::rxIs(.cand, "rxSymInvCholEnv")) {
        .omChk <- tryCatch(as.matrix(.cand$omega), error = function(.) NULL)
        if (!is.null(.omChk) && identical(dim(.omChk), dim(as.matrix(Om))) &&
              isTRUE(all.equal(unname(.omChk), unname(as.matrix(Om)), tolerance = 1e-10))) {
          .rxInv <- .cand
        }
      }
    }
    if (is.null(.rxInv)) .rxInv <- rxode2::rxSymInvCholCreate(mat = Om, diag.xform = .diagXform)
    # `$.rxSymInvCholEnv` dispatches to the C rxSymInvCholEnvCalculate; d.omegaInv
    # is the list of dOmega^-1/d(chol theta_k), tr.28 = 0.5*tr(dOmega^-1_k Omega).
    .dOi <- .rxInv$d.omegaInv
    .tr28 <- .rxInv$tr.28
    if (is.null(.dOi) || is.null(.tr28)) return(NULL)
    list(dOi = .dOi, tr28 = as.numeric(.tr28),
         names = paste0("om.chol.", seq_along(.dOi)))
  }, error = function(e) NULL)
}

#' Theta names excluded from the outer optimizer's free-parameter set by the
#' mu-referenced (lin/irls) regression -- mirrors inner.cpp isMuGroupSkip: the
#' mu-group thetas plus every mu-group covariate coefficient (bounded ones are
#' regression-updated with clamping too).  Index arrays are 0-based (see
#' `.muRefCppGroupSetup`).
#' @noRd
.foceiMuSkipThetaNames <- function(ui, thNames) {
  if (identical(rxode2::rxGetControl(ui, "muModel", "none"), "none")) {
    return(character(0))
  }
  .g <- as.integer(rxode2::rxGetControl(ui, "foceiMuGroupTheta", integer(0)))
  .ct <- as.integer(rxode2::rxGetControl(ui, "foceiMuGroupCovTheta", integer(0)))
  thNames[c(.g, .ct) + 1L]
}

