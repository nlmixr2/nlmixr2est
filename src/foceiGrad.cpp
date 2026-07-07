// C++/Armadillo port of the FOCEI per-subject analytic outer-gradient assembly
// (R oracle: .foceiAnalyticSubjectGrad in R/foceiGradAnalytic.R).  The R driver
// solves the augmented sensitivity model and evaluates the per-observation error
// coefficients (r1/r2/p/p1 and, per residual sigma, rf/ps/rs); this kernel does
// the O(neta^2 * nobs) tensor contractions (H/Ht/N/dHtD/etaP) and returns the
// length-np gradient of the objective (-2*logLik) plus etaP (d eta*/d p, Eq 46).
//
// Param order matches the R: nth structural theta, nsg residual sigma, nom Omega
// (Cholesky) params.  dirTh is 1-based (a mu-ref theta reuses its eta's direction).
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "rxomp.h"
using namespace arma;

// [[Rcpp::export]]
Rcpp::List foceiSubjectGradFocei_(const arma::mat& a,       // nobs x ndir  (d f / d dir)
                                  const arma::cube& A,       // nobs x ndir x ndir (2nd order)
                                  const arma::vec& r1, const arma::vec& r2,  // nobs (rho f-derivs)
                                  const arma::vec& p,  const arma::vec& p1,   // nobs (det p, dp/df)
                                  const arma::mat& perRf,    // nobs x nsg  (d2 rho / df dsig)
                                  const arma::mat& perPs,    // nobs x nsg  (d p / dsig)
                                  const arma::mat& perRs,    // nobs x nsg  (d rho / dsig)
                                  const arma::vec& ehat,     // neta (EBE)
                                  const arma::mat& Oi,       // neta x neta (Omega^-1)
                                  const arma::cube& dOiEst,  // neta x neta x nom (est-scale dOmega^-1)
                                  const arma::vec& tr28,     // nom (0.5 tr(dOmega^-1 Omega))
                                  int neta, int nth, int nsg, int nom,
                                  const arma::ivec& dirTh) { // nth (1-based direction per theta)
  const int nobs = (int)a.n_rows;
  const int ndir = (int)a.n_cols;
  const int np = nth + nsg + nom;

  // Inner Hessian H = d2l/deta2 (for etaP) and Laplace determinant Hessian Ht.
  mat H = Oi, Ht = Oi;
  for (int l = 0; l < neta; l++) {
    for (int m = 0; m < neta; m++) {
      double sh = 0.0, sht = 0.0;
      for (int o = 0; o < nobs; o++) {
        sh  += r2[o] * a(o, l) * a(o, m) + r1[o] * A(o, l, m);
        sht += p[o]  * a(o, l) * a(o, m);
      }
      H(l, m) += sh; Ht(l, m) += sht;
    }
  }
  mat HiM = inv(H), Hti = inv(Ht);

  // N[l,d] = d2l/(deta_l ddir_d)
  mat N(neta, ndir, fill::zeros);
  for (int l = 0; l < neta; l++)
    for (int d = 0; d < ndir; d++) {
      double s = 0.0;
      for (int o = 0; o < nobs; o++) s += r2[o] * a(o, l) * a(o, d) + r1[o] * A(o, l, d);
      N(l, d) = s;
    }

  // dHtD[s] = dHt/d(direction s); eta-directions (s < neta) are the moving mode.
  std::vector<mat> dHtD(ndir);
  for (int s = 0; s < ndir; s++) {
    mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++)
      for (int m = 0; m < neta; m++) {
        double v = 0.0;
        for (int o = 0; o < nobs; o++)
          v += p1[o] * a(o, s) * a(o, l) * a(o, m) + p[o] * A(o, l, s) * a(o, m) + p[o] * a(o, l) * A(o, m, s);
        D(l, m) = v;
      }
    dHtD[s] = D;
  }

  // ouAA(v)[l,m] = sum_o v_o a_l a_m
  auto ouAA = [&](const vec& v) {
    mat M(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++)
      for (int m = 0; m < neta; m++) {
        double s = 0.0;
        for (int o = 0; o < nobs; o++) s += v[o] * a(o, l) * a(o, m);
        M(l, m) = s;
      }
    return M;
  };

  // M_p = d2l/(deta dp): th -> N[,dir]; sg -> a'(d2rho/df dsig); om -> dOmega^-1 eta
  auto Mcol = [&](int pp) {
    vec r(neta, fill::zeros);
    if (pp < nth) {
      r = N.col(dirTh[pp] - 1);
    } else if (pp < nth + nsg) {
      int j = pp - nth;
      for (int l = 0; l < neta; l++) {
        double s = 0.0;
        for (int o = 0; o < nobs; o++) s += a(o, l) * perRf(o, j);
        r[l] = s;
      }
    } else {
      r = dOiEst.slice(pp - nth - nsg) * ehat;
    }
    return r;
  };

  mat etaP(neta, np, fill::zeros);
  for (int pp = 0; pp < np; pp++) etaP.col(pp) = -HiM * Mcol(pp);

  vec g(np, fill::zeros);
  for (int pp = 0; pp < np; pp++) {
    mat dHtStar;
    double dPhi;
    if (pp < nth) {
      dHtStar = dHtD[dirTh[pp] - 1];
      double s = 0.0;
      for (int o = 0; o < nobs; o++) s += r1[o] * a(o, dirTh[pp] - 1);
      dPhi = s;
    } else if (pp < nth + nsg) {
      int j = pp - nth;
      dHtStar = ouAA(perPs.col(j));
      double s = 0.0;
      for (int o = 0; o < nobs; o++) s += perRs(o, j);
      dPhi = s;
    } else {
      int k = pp - nth - nsg;
      dHtStar = dOiEst.slice(k);
      dPhi = 0.5 * as_scalar(ehat.t() * dOiEst.slice(k) * ehat) - tr28[k];
    }
    for (int l = 0; l < neta; l++) dHtStar += etaP(l, pp) * dHtD[l];
    g[pp] = 2.0 * dPhi + trace(Hti * dHtStar);
  }

  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaP);
}

// (f,R) FOCEI per-subject outer gradient (R oracle: .foceiAnalyticSubjectGradFR).
// The prediction f and the variance R are independent solved quantities: a/A are the
// prediction sensitivities, aR/AR the variance sensitivities, and a residual sigma is a
// pseudo-direction (df/dsigma=0, dR/dsigma=Rsig).  The rho(f,R,y) partials are
// model-independent closed forms computed here from f/y/R, so ANY variance structure works.
// Shared per-subject core, called both from the single-subject export (oracle) and the
// batched OpenMP driver foceiGradAllFR_; writes g_out (np) and etaP_out (neta x np).
static void foceiGradSubjectFR_(const arma::mat& a, const arma::cube& A,
                                const arma::mat& aR, const arma::cube& AR,
                                const arma::mat& Rsig, const arma::cube& RsigDir,
                                const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                                const arma::vec& ehat, const arma::mat& Oi,
                                const arma::cube& dOiEst, const arma::vec& tr28,
                                int neta, int nth, int nsg, int nom,
                                const arma::ivec& dirTh, const arma::ivec& sigCol,
                                arma::vec& g_out, arma::mat& etaP_out) {
  const int nobs = (int)a.n_rows;
  const int ndir = (int)a.n_cols;
  const int np = nth + nsg + nom;
  // per-observation rho(f,R,y) partials
  vec res = yv - fv;
  vec rf = -res / Rv, rR = 0.5 * (1.0 / Rv - square(res) / square(Rv));
  vec rff = 1.0 / Rv, rfR = res / square(Rv), rRR = 0.5 * (-1.0 / square(Rv) + 2.0 * square(res) / pow(Rv, 3));
  vec eff = 1.0 / Rv, eRR = 0.5 / square(Rv);
  // exact inner Hessian H (eta x eta), N (eta x dir), determinant Ht
  mat H = Oi, Ht = Oi, N(neta, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) {
    for (int m = 0; m < neta; m++) {
      double sh = 0.0, sht = 0.0;
      for (int o = 0; o < nobs; o++) {
        sh += rff[o] * a(o, l) * a(o, m) + rfR[o] * (a(o, l) * aR(o, m) + aR(o, l) * a(o, m)) +
          rRR[o] * aR(o, l) * aR(o, m) + rf[o] * A(o, l, m) + rR[o] * AR(o, l, m);
        sht += eff[o] * a(o, l) * a(o, m) + eRR[o] * aR(o, l) * aR(o, m);
      }
      H(l, m) += sh; Ht(l, m) += sht;
    }
    for (int d = 0; d < ndir; d++) {
      double s = 0.0;
      for (int o = 0; o < nobs; o++)
        s += rff[o] * a(o, l) * a(o, d) + rfR[o] * (a(o, l) * aR(o, d) + aR(o, l) * a(o, d)) +
          rRR[o] * aR(o, l) * aR(o, d) + rf[o] * A(o, l, d) + rR[o] * AR(o, l, d);
      N(l, d) = s;
    }
  }
  mat HiM = inv(H), Hti = inv(Ht);
  // eta x sigma block (df/dsigma=0)
  mat Nsg(neta, nsg, fill::zeros);
  for (int l = 0; l < neta; l++)
    for (int k = 0; k < nsg; k++) {
      int c = sigCol[k] - 1; double s = 0.0;
      for (int o = 0; o < nobs; o++)
        s += rfR[o] * a(o, l) * Rsig(o, c) + rRR[o] * aR(o, l) * Rsig(o, c) + rR[o] * RsigDir(o, l, c);
      Nsg(l, k) = s;
    }
  // dHt/d(direction s) (moving mode s=eta and explicit theta columns)
  std::vector<mat> dHtD(ndir);
  for (int s = 0; s < ndir; s++) {
    mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++)
      for (int m = 0; m < neta; m++) {
        double v = 0.0;
        for (int o = 0; o < nobs; o++)
          v += -aR(o, s) / (Rv[o] * Rv[o]) * a(o, l) * a(o, m) + eff[o] * (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) +
            -aR(o, s) / pow(Rv[o], 3) * aR(o, l) * aR(o, m) + eRR[o] * (AR(o, l, s) * aR(o, m) + aR(o, l) * AR(o, m, s));
        D(l, m) = v;
      }
    dHtD[s] = D;
  }
  // dHt/dsigma_k (df/dsigma=0)
  std::vector<mat> dHtSg(nsg);
  for (int k = 0; k < nsg; k++) {
    int c = sigCol[k] - 1; mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++)
      for (int m = 0; m < neta; m++) {
        double v = 0.0;
        for (int o = 0; o < nobs; o++)
          v += -Rsig(o, c) / (Rv[o] * Rv[o]) * a(o, l) * a(o, m) - Rsig(o, c) / pow(Rv[o], 3) * aR(o, l) * aR(o, m) +
            eRR[o] * (RsigDir(o, l, c) * aR(o, m) + aR(o, l) * RsigDir(o, m, c));
        D(l, m) = v;
      }
    dHtSg[k] = D;
  }
  auto Mcol = [&](int pp) {
    vec r(neta, fill::zeros);
    if (pp < nth) r = N.col(dirTh[pp] - 1);
    else if (pp < nth + nsg) r = Nsg.col(pp - nth);
    else r = dOiEst.slice(pp - nth - nsg) * ehat;
    return r;
  };
  mat etaP(neta, np, fill::zeros);
  for (int pp = 0; pp < np; pp++) etaP.col(pp) = -HiM * Mcol(pp);
  vec g(np, fill::zeros);
  for (int pp = 0; pp < np; pp++) {
    mat dHtStar; double dPhi;
    if (pp < nth) {
      int d = dirTh[pp] - 1; dHtStar = dHtD[d];
      double s = 0.0; for (int o = 0; o < nobs; o++) s += rf[o] * a(o, d) + rR[o] * aR(o, d);
      dPhi = s;
    } else if (pp < nth + nsg) {
      int c = sigCol[pp - nth] - 1; dHtStar = dHtSg[pp - nth];
      double s = 0.0; for (int o = 0; o < nobs; o++) s += rR[o] * Rsig(o, c);
      dPhi = s;
    } else {
      int k = pp - nth - nsg; dHtStar = dOiEst.slice(k);
      dPhi = 0.5 * as_scalar(ehat.t() * dOiEst.slice(k) * ehat) - tr28[k];
    }
    for (int l = 0; l < neta; l++) dHtStar += etaP(l, pp) * dHtD[l];
    g[pp] = 2.0 * dPhi + trace(Hti * dHtStar);
  }
  g_out = g; etaP_out = etaP;
}

// Single-subject export (oracle / R fallback): thin wrapper over foceiGradSubjectFR_.
// [[Rcpp::export]]
Rcpp::List foceiSubjectGradFR_(const arma::mat& a, const arma::cube& A,
                               const arma::mat& aR, const arma::cube& AR,
                               const arma::mat& Rsig, const arma::cube& RsigDir,
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                               const arma::vec& ehat, const arma::mat& Oi,
                               const arma::cube& dOiEst, const arma::vec& tr28,
                               int neta, int nth, int nsg, int nom,
                               const arma::ivec& dirTh, const arma::ivec& sigCol) {
  vec g; mat etaP;
  foceiGradSubjectFR_(a, A, aR, AR, Rsig, RsigDir, fv, yv, Rv, ehat, Oi, dOiEst, tr28,
                      neta, nth, nsg, nom, dirTh, sigCol, g, etaP);
  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaP);
}

// Batched (f,R) FOCEI outer gradient over ALL subjects in one OpenMP-parallel C++ call:
// removes the per-subject R<->C++ round-trip.  Sensitivities are concatenated over
// observations (obsOffset[i]..obsOffset[i+1]-1 are subject i's rows); ehat is nsub x neta.
// Returns the summed gradient g (np) and the per-subject etaP as a cube (neta x np x nsub).
// [[Rcpp::export]]
Rcpp::List foceiGradAllFR_(const arma::mat& a, const arma::cube& A,
                           const arma::mat& aR, const arma::cube& AR,
                           const arma::mat& Rsig, const arma::cube& RsigDir,
                           const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                           const arma::mat& ehat, const arma::ivec& obsOffset,
                           const arma::mat& Oi, const arma::cube& dOiEst, const arma::vec& tr28,
                           int neta, int nth, int nsg, int nom,
                           const arma::ivec& dirTh, const arma::ivec& sigCol, int ncores) {
  const int nsub = (int)ehat.n_rows;
  const int np = nth + nsg + nom;
  const int ndir = (int)a.n_cols;
  mat gmat(np, nsub, fill::zeros);
  cube etaPall(neta, np, nsub, fill::zeros);
  const bool hasSig = (Rsig.n_cols > 0);
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < nsub; i++) {
    int o0 = obsOffset[i], o1 = obsOffset[i + 1] - 1;
    mat ai = a.rows(o0, o1), aRi = aR.rows(o0, o1);
    cube Ai = A.rows(o0, o1), ARi = AR.rows(o0, o1);
    mat Rsigi = hasSig ? mat(Rsig.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
    cube RsigDiri = hasSig ? cube(RsigDir.rows(o0, o1)) : cube(o1 - o0 + 1, ndir, 0);
    vec gi; mat etaPi;
    foceiGradSubjectFR_(ai, Ai, aRi, ARi, Rsigi, RsigDiri, fv.subvec(o0, o1), yv.subvec(o0, o1),
                        Rv.subvec(o0, o1), ehat.row(i).t(), Oi, dOiEst, tr28,
                        neta, nth, nsg, nom, dirTh, sigCol, gi, etaPi);
    gmat.col(i) = gi; etaPall.slice(i) = etaPi;
  }
  vec g = sum(gmat, 1);
  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaPall);
}

// (f,R) FOCE per-subject outer gradient (oracle: .foceiAnalyticSubjectGradFoceFR).  The
// inner problem is interaction-free (q0=-(y-f)/R0, q1=1/R0) with a frozen variance R0: the
// eta-block (gPhi/Hf/Ht/dHtD) uses aRe (0 for nonmem, live E$aR for foce+) and the parameter
// columns use aRc (E0's dR0/ddir for nonmem, live E$aR for foce+).  The determinant is the
// Gauss-Newton Ht=Omega^-1+sum(a a/R0).  `fp` = foce+ (1) vs nonmem (0): nonmem adds the
// aRc a0-chain to dHt/dtheta (dHtD used aRe=0 there).  No 3rd-order tensor (gradient only).
// [[Rcpp::export]]
Rcpp::List foceiSubjectGradFoceFR_(const arma::mat& a, const arma::cube& A,
                                   const arma::mat& aRe, const arma::mat& aRc,
                                   const arma::mat& R0sig,
                                   const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                                   const arma::vec& ehat, const arma::mat& Oi,
                                   const arma::cube& dOiEst, const arma::vec& tr28,
                                   int neta, int nth, int nsg, int nom,
                                   const arma::ivec& dirTh, const arma::ivec& sigCol, int fp) {
  const int nobs = (int)a.n_rows;
  const int ndir = (int)a.n_cols;
  const int np = nth + nsg + nom;
  vec res = yv - fv;
  vec rho_f = -res / R0v, rho_R = 0.5 * (1.0 / R0v - square(res) / square(R0v));
  vec q0 = rho_f, q1 = 1.0 / R0v;
  vec iR = 1.0 / R0v, iR2 = square(iR);
  vec gPhi = Oi * ehat;
  for (int l = 0; l < neta; l++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += rho_f[o] * a(o, l) + rho_R[o] * aRe(o, l); gPhi[l] += s; }
  mat Hf = Oi, Ht = Oi, Nf(neta, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) {
    for (int m = 0; m < neta; m++) { double sh = 0.0, st = 0.0;
      for (int o = 0; o < nobs; o++) { sh += q1[o] * a(o, l) * a(o, m) + q0[o] * A(o, l, m); st += a(o, l) * a(o, m) * iR[o]; }
      Hf(l, m) += sh; Ht(l, m) += st; }
    for (int d = 0; d < ndir; d++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += q1[o] * a(o, l) * a(o, d) + q0[o] * A(o, l, d); Nf(l, d) = s; }
  }
  mat HfInv = inv(Hf), Hti = inv(Ht);
  std::vector<mat> dHtD(ndir);
  for (int s = 0; s < ndir; s++) { mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
      for (int o = 0; o < nobs; o++) v += -aRe(o, s) * iR2[o] * a(o, l) * a(o, m) + (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) * iR[o];
      D(l, m) = v; }
    dHtD[s] = D; }
  vec Cen(neta); for (int l = 0; l < neta; l++) Cen[l] = 0.5 * trace(Hti * dHtD[l]);
  auto ouRc = [&](const arma::vec& v) { mat M(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += v[o] * a(o, l) * a(o, m); M(l, m) = s; } return M; };
  auto typ = [&](int p) { return p < nth ? 0 : (p < nth + nsg ? 1 : 2); };
  auto McolEBE = [&](int p) -> vec { int t = typ(p);
    if (t == 0) { int d = dirTh[p] - 1; vec r = Nf.col(d);
      for (int l = 0; l < neta; l++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += a(o, l) * (res[o] * iR2[o]) * aRc(o, d); r[l] += s; } return r; }
    if (t == 1) { int c = sigCol[p - nth] - 1; vec r(neta);
      for (int l = 0; l < neta; l++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += a(o, l) * (res[o] * iR2[o]) * R0sig(o, c); r[l] = s; } return r; }
    return vec(dOiEst.slice(p - nth - nsg) * ehat); };
  auto dHt_p = [&](int p) -> mat { int t = typ(p);
    if (t == 0) { int d = dirTh[p] - 1; mat D = dHtD[d]; if (!fp) D += ouRc(-aRc.col(d) % iR2); return D; }
    if (t == 1) { int c = sigCol[p - nth] - 1; return ouRc(-R0sig.col(c) % iR2); }
    return dOiEst.slice(p - nth - nsg); };
  auto dPhiExplicit = [&](int p) -> double { int t = typ(p);
    if (t == 0) { int d = dirTh[p] - 1; double s = 0.0; for (int o = 0; o < nobs; o++) s += rho_f[o] * a(o, d) + rho_R[o] * aRc(o, d); return s; }
    if (t == 1) { int c = sigCol[p - nth] - 1; double s = 0.0; for (int o = 0; o < nobs; o++) s += rho_R[o] * R0sig(o, c); return s; }
    int k = p - nth - nsg; return 0.5 * as_scalar(ehat.t() * dOiEst.slice(k) * ehat) - tr28[k]; };
  mat etaP(neta, np, fill::zeros);
  for (int p = 0; p < np; p++) etaP.col(p) = -HfInv * McolEBE(p);
  vec g(np, fill::zeros);
  for (int p = 0; p < np; p++) {
    double v = dPhiExplicit(p) + 0.5 * trace(Hti * dHt_p(p));
    for (int l = 0; l < neta; l++) v += (gPhi[l] + Cen[l]) * etaP(l, p);
    g[p] = 2.0 * v;
  }
  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaP);
}

// (f,R) FOCEI per-subject observed-information R (oracle: .foceiAnalyticSubjectRFR).
// Analytic 1st/2nd-order sensitivities a/A (prediction) + aR/AR (variance); the 3rd-order
// tensors Ath/AthR come from Shi-FD and are passed reshaped to cubes (nobs, ndir, ndir*ndir):
// Ath[o,l,s,t] == Ath(o, l, s + t*ndir).  Every non-Omega param is a direction (dirP,
// 1-based); a residual sigma direction has a=A=Ath=0.  Omega derivatives: dOi
// (neta,neta,nom), d2Oi (neta,neta,nom*nom) with slice a*nom+b, d2LD (nom,nom).
// [[Rcpp::export]]
arma::mat foceiSubjectRFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                           const arma::mat& aR, const arma::cube& AR, const arma::cube& AthR,
                           const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                           const arma::vec& ehat, const arma::mat& Oi,
                           const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                           int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP) {
  const int nobs = (int)a.n_rows;
  const int np = ndirP + nom;
  vec res = yv - fv;
  vec rf = -res / Rv, rR = 0.5 * (1.0 / Rv - square(res) / square(Rv));
  vec rff = 1.0 / Rv, rfR = res / square(Rv), rRR = 0.5 * (-1.0 / square(Rv) + 2.0 * square(res) / pow(Rv, 3));
  vec rffR = -1.0 / square(Rv), rfRR = -2.0 * res / pow(Rv, 3), rRRR = 0.5 * (2.0 / pow(Rv, 3) - 6.0 * square(res) / pow(Rv, 4));
  vec iR = 1.0 / Rv, iR2 = square(iR), iR3 = pow(iR, 3), iR4 = pow(iR, 4);
  auto Ai = [&](const arma::cube& T, int o, int l, int s, int t) { return T(o, l, s + t * ndir); };
  // Gdd: (f,R) 2nd total derivative of the density between two directions
  auto Gdd = [&](int da, int db) {
    double s = 0.0;
    for (int o = 0; o < nobs; o++)
      s += rff[o] * a(o, da) * a(o, db) + rfR[o] * (a(o, da) * aR(o, db) + aR(o, da) * a(o, db)) +
        rRR[o] * aR(o, da) * aR(o, db) + rf[o] * A(o, da, db) + rR[o] * AR(o, da, db);
    return s;
  };
  mat H = Oi; for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) H(l, m) += Gdd(l, m);
  mat HiM = inv(H);
  mat N(neta, ndir, fill::zeros); for (int l = 0; l < neta; l++) for (int d = 0; d < ndir; d++) N(l, d) = Gdd(l, d);
  // Tn[l,s,t] = d2(rf a_l + rR aR_l)/ddir_s ddir_t  -> cube(neta, ndir, ndir), Tn(l,s,t)=Tn.slice(t)(l,s)
  cube Tn(neta, ndir, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) for (int s = 0; s < ndir; s++) for (int t = 0; t < ndir; t++) {
    double v = 0.0;
    for (int o = 0; o < nobs; o++) {
      double us = rff[o] * a(o, s) + rfR[o] * aR(o, s), ut = rff[o] * a(o, t) + rfR[o] * aR(o, t);
      double ust = rffR[o] * (a(o, s) * aR(o, t) + aR(o, s) * a(o, t)) + rfRR[o] * aR(o, s) * aR(o, t) +
        rff[o] * A(o, s, t) + rfR[o] * AR(o, s, t);
      double ws = rfR[o] * a(o, s) + rRR[o] * aR(o, s), wt = rfR[o] * a(o, t) + rRR[o] * aR(o, t);
      double wst = rffR[o] * a(o, s) * a(o, t) + rfRR[o] * (a(o, s) * aR(o, t) + aR(o, s) * a(o, t)) +
        rRRR[o] * aR(o, s) * aR(o, t) + rfR[o] * A(o, s, t) + rRR[o] * AR(o, s, t);
      v += ust * a(o, l) + us * A(o, l, t) + ut * A(o, l, s) + rf[o] * Ai(Ath, o, l, s, t) +
        wst * aR(o, l) + ws * AR(o, l, t) + wt * AR(o, l, s) + rR[o] * Ai(AthR, o, l, s, t);
    }
    Tn(l, s, t) = v;
  }
  // determinant Ht = Oi + sum(a a / R + 0.5 aR aR / R^2); dHtD / d2HtDD by (f,R) chain
  mat Ht = Oi;
  for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) {
    double v = 0.0; for (int o = 0; o < nobs; o++) v += a(o, l) * a(o, m) * iR[o] + 0.5 * aR(o, l) * aR(o, m) * iR2[o];
    Ht(l, m) += v;
  }
  mat Hti = inv(Ht);
  std::vector<mat> dHtD(ndir);
  for (int s = 0; s < ndir; s++) { mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
      for (int o = 0; o < nobs; o++)
        v += (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) * iR[o] - a(o, l) * a(o, m) * aR(o, s) * iR2[o] +
          0.5 * (AR(o, l, s) * aR(o, m) + aR(o, l) * AR(o, m, s)) * iR2[o] - aR(o, l) * aR(o, m) * aR(o, s) * iR3[o];
      D(l, m) = v; }
    dHtD[s] = D; }
  std::vector<std::vector<mat> > d2HtDD(ndir, std::vector<mat>(ndir));
  for (int s = 0; s < ndir; s++) for (int t = 0; t < ndir; t++) { mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
      for (int o = 0; o < nobs; o++)
        v += (Ai(Ath, o, l, s, t) * a(o, m) + A(o, l, s) * A(o, m, t) + A(o, l, t) * A(o, m, s) + a(o, l) * Ai(Ath, o, m, s, t)) * iR[o] -
          (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) * aR(o, t) * iR2[o] -
          (A(o, l, t) * a(o, m) + a(o, l) * A(o, m, t)) * aR(o, s) * iR2[o] - a(o, l) * a(o, m) * AR(o, s, t) * iR2[o] +
          2.0 * a(o, l) * a(o, m) * aR(o, s) * aR(o, t) * iR3[o] +
          0.5 * (Ai(AthR, o, l, s, t) * aR(o, m) + AR(o, l, s) * AR(o, m, t) + AR(o, l, t) * AR(o, m, s) + aR(o, l) * Ai(AthR, o, m, s, t)) * iR2[o] -
          (AR(o, l, s) * aR(o, m) + aR(o, l) * AR(o, m, s)) * aR(o, t) * iR3[o] -
          (AR(o, l, t) * aR(o, m) + aR(o, l) * AR(o, m, t)) * aR(o, s) * iR3[o] - aR(o, l) * aR(o, m) * AR(o, s, t) * iR3[o] +
          3.0 * aR(o, l) * aR(o, m) * aR(o, s) * aR(o, t) * iR4[o];
      D(l, m) = v; }
    d2HtDD[s][t] = D; }
  vec Cen(neta); for (int l = 0; l < neta; l++) Cen[l] = 0.5 * trace(Hti * dHtD[l]);
  mat Cee(neta, neta, fill::zeros);
  for (int s = 0; s < neta; s++) for (int t = 0; t < neta; t++)
    Cee(s, t) = 0.5 * (trace(Hti * d2HtDD[s][t]) - trace(Hti * dHtD[s] * Hti * dHtD[t]));
  auto isDir = [&](int p) { return p < ndirP; };
  auto dOf = [&](int p) { return dirP[p] - 1; };
  auto omc = [&](int p) { return p - ndirP; };
  auto Mcol = [&](int p) { vec r(neta);
    if (isDir(p)) r = N.col(dOf(p)); else r = dOi.slice(omc(p)) * ehat; return r; };
  auto dHt_p = [&](int p) -> mat { if (isDir(p)) return dHtD[dOf(p)]; return dOi.slice(omc(p)); };
  auto d2HtEtaP = [&](int p, int l) -> mat { if (isDir(p)) return d2HtDD[dOf(p)][l]; return zeros<mat>(neta, neta); };
  auto d2Ht_pp = [&](int aa, int bb) -> mat {
    if (isDir(aa) && isDir(bb)) return d2HtDD[dOf(aa)][dOf(bb)];
    if (!isDir(aa) && !isDir(bb)) return d2Oi.slice(omc(aa) * nom + omc(bb));
    return zeros<mat>(neta, neta); };
  auto Smat = [&](int p) { mat M(neta, ndir, fill::zeros);
    if (!isDir(p)) { M.cols(0, neta - 1) = dOi.slice(omc(p)); return M; }
    for (int l = 0; l < neta; l++) for (int s = 0; s < ndir; s++) M(l, s) = Tn(l, dOf(p), s); return M; };
  auto Svec = [&](int aa, int bb) { vec v(neta, fill::zeros);
    if (isDir(aa) && isDir(bb)) { for (int l = 0; l < neta; l++) v[l] = Tn(l, dOf(aa), dOf(bb)); return v; }
    if (!isDir(aa) && !isDir(bb)) { v = d2Oi.slice(omc(aa) * nom + omc(bb)) * ehat; return v; }
    return v; };
  auto d2Phi = [&](int aa, int bb) {
    if (isDir(aa) && isDir(bb)) return Gdd(dOf(aa), dOf(bb));
    if (!isDir(aa) && !isDir(bb)) return 0.5 * as_scalar(ehat.t() * d2Oi.slice(omc(aa) * nom + omc(bb)) * ehat) + 0.5 * d2LD(omc(aa), omc(bb));
    return 0.0; };
  std::vector<vec> Mcols(np); for (int p = 0; p < np; p++) Mcols[p] = Mcol(p);
  mat etaP(neta, np); for (int p = 0; p < np; p++) etaP.col(p) = -HiM * Mcols[p];
  auto eta2 = [&](int aa, int bb) {
    mat SmA = Smat(aa).cols(0, neta - 1), SmB = Smat(bb).cols(0, neta - 1);
    vec b = Svec(aa, bb) + SmA * etaP.col(bb) + SmB * etaP.col(aa);
    for (int l = 0; l < neta; l++) { double acc = 0.0;
      for (int s = 0; s < neta; s++) for (int t = 0; t < neta; t++) acc += etaP(s, aa) * Tn(l, s, t) * etaP(t, bb);
      b[l] += acc; }
    return vec(-HiM * b); };
  auto Cpe = [&](int p, int l) { return 0.5 * (trace(Hti * d2HtEtaP(p, l)) - trace(Hti * dHt_p(p) * Hti * dHtD[l])); };
  auto Cpp = [&](int aa, int bb) { return 0.5 * (trace(Hti * d2Ht_pp(aa, bb)) - trace(Hti * dHt_p(aa) * Hti * dHt_p(bb))); };
  std::vector<vec> CpeRow(np);
  for (int p = 0; p < np; p++) { vec r(neta); for (int l = 0; l < neta; l++) r[l] = Cpe(p, l); CpeRow[p] = r; }
  mat R(np, np, fill::zeros);
  for (int aa = 0; aa < np; aa++) for (int bb = aa; bb < np; bb++) {
    double dat = d2Phi(aa, bb) - as_scalar(Mcols[aa].t() * HiM * Mcols[bb]);
    vec e2 = eta2(aa, bb);
    double ld = Cpp(aa, bb) + dot(CpeRow[aa], etaP.col(bb)) + dot(CpeRow[bb], etaP.col(aa)) +
      as_scalar(etaP.col(aa).t() * Cee * etaP.col(bb)) + dot(Cen, e2);
    R(aa, bb) = R(bb, aa) = dat + ld;
  }
  return R;
}

// (f,R) FOCE per-subject observed-information R (oracle: .foceiAnalyticSubjectRfoceFR).
// Interaction-free inner (Hf = Oi + sum(q1 a a + q0 A), q0 = -(y-f)/R0, q1 = 1/R0) with a
// frozen variance R0 and the non-envelope assembly (Phi_eta = S_FOCE ~ 0 at the EBE but the
// log-determinant's eta-gradient is not).  R0's theta-chain enters the parameter columns via
// aRc/ARc while the eta-block stays frozen: aRe/ARe drive the eta-block (0 for nonmem, the
// live E$aR/E$AR for foce+) and aRc/ARc the parameter columns (E0's dR0/ddir, d2R0/ddir2 for
// nonmem, the same live E for foce+).  Ath is reshaped as in foceiSubjectRFR_; a sigma
// direction has a=A=Ath=0 (only aRc/ARc).
// [[Rcpp::export]]
arma::mat foceiSubjectRfoceFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                               const arma::mat& aRe, const arma::mat& aRc,
                               const arma::cube& ARe, const arma::cube& ARc,
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                               const arma::vec& ehat, const arma::mat& Oi,
                               const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                               int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP) {
  const int nobs = (int)a.n_rows;
  const int np = ndirP + nom;
  vec res = yv - fv;
  vec rf = -res / R0v, rR = 0.5 * (1.0 / R0v - square(res) / square(R0v));
  vec rff = 1.0 / R0v, rfR = res / square(R0v), rRR = 0.5 * (-1.0 / square(R0v) + 2.0 * square(res) / pow(R0v, 3));
  vec q0 = rf, q1 = rff;                              // interaction-free inner
  vec iR = 1.0 / R0v, iR2 = square(iR), iR3 = pow(iR, 3);
  auto Ai = [&](const arma::cube& T, int o, int l, int s, int t) { return T(o, l, s + t * ndir); };
  // ---- Phi (data) tensors: H = Phi_etaeta, gPhi = Phi_eta (aRe eta-block) ----
  mat H = Oi;
  for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
    for (int o = 0; o < nobs; o++)
      v += rff[o] * a(o, l) * a(o, m) + rfR[o] * (a(o, l) * aRe(o, m) + aRe(o, l) * a(o, m)) +
        rRR[o] * aRe(o, l) * aRe(o, m) + rf[o] * A(o, l, m) + rR[o] * ARe(o, l, m);
    H(l, m) += v; }
  vec gPhi = Oi * ehat;
  for (int l = 0; l < neta; l++) { double v = 0.0; for (int o = 0; o < nobs; o++) v += rf[o] * a(o, l) + rR[o] * aRe(o, l); gPhi[l] += v; }
  // ---- FOCE inner (EBE) tensors: interaction-free q-based Hf/Nf/Tnf ----
  mat Hf = Oi; mat Nf(neta, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) {
    for (int m = 0; m < neta; m++) { double v = 0.0; for (int o = 0; o < nobs; o++) v += q1[o] * a(o, l) * a(o, m) + q0[o] * A(o, l, m); Hf(l, m) += v; }
    for (int d = 0; d < ndir; d++) { double v = 0.0; for (int o = 0; o < nobs; o++) v += q1[o] * a(o, l) * a(o, d) + q0[o] * A(o, l, d); Nf(l, d) = v; } }
  mat HfInv = inv(Hf);
  cube Tnf(neta, ndir, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) for (int s = 0; s < ndir; s++) for (int t = 0; t < ndir; t++) { double v = 0.0;
    for (int o = 0; o < nobs; o++) v += q1[o] * (A(o, l, s) * a(o, t) + A(o, l, t) * a(o, s) + A(o, s, t) * a(o, l)) + q0[o] * Ai(Ath, o, l, s, t);
    Tnf(l, s, t) = v; }
  // ---- determinant Ht = Oi + sum(a a / R0) (interaction-free) + its derivatives ----
  mat Ht = Oi; for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
    for (int o = 0; o < nobs; o++) v += a(o, l) * a(o, m) * iR[o]; Ht(l, m) += v; }
  mat Hti = inv(Ht);
  auto dHtDir = [&](int s, const arma::mat& aRv) -> mat { mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
      for (int o = 0; o < nobs; o++) v += (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) * iR[o] - a(o, l) * a(o, m) * aRv(o, s) * iR2[o];
      D(l, m) = v; } return D; };
  auto d2HtDir = [&](int s, int t, const arma::mat& aRvS, const arma::mat& aRvT, const arma::cube& ARv) -> mat {
    mat D(neta, neta, fill::zeros);
    for (int l = 0; l < neta; l++) for (int m = 0; m < neta; m++) { double v = 0.0;
      for (int o = 0; o < nobs; o++) v += (Ai(Ath, o, l, s, t) * a(o, m) + A(o, l, s) * A(o, m, t) + A(o, l, t) * A(o, m, s) + a(o, l) * Ai(Ath, o, m, s, t)) * iR[o] -
        (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) * aRvT(o, t) * iR2[o] -
        (A(o, l, t) * a(o, m) + a(o, l) * A(o, m, t)) * aRvS(o, s) * iR2[o] - a(o, l) * a(o, m) * ARv(o, s, t) * iR2[o] +
        2.0 * a(o, l) * a(o, m) * aRvS(o, s) * aRvT(o, t) * iR3[o];
      D(l, m) = v; } return D; };
  std::vector<mat> dHtE(neta); for (int l = 0; l < neta; l++) dHtE[l] = dHtDir(l, aRe);
  std::vector<std::vector<mat> > d2HtEE(neta, std::vector<mat>(neta));
  for (int s = 0; s < neta; s++) for (int t = 0; t < neta; t++) d2HtEE[s][t] = d2HtDir(s, t, aRe, aRe, ARe);
  vec Cen(neta); for (int l = 0; l < neta; l++) Cen[l] = 0.5 * trace(Hti * dHtE[l]);
  mat Cee(neta, neta, fill::zeros);
  for (int s = 0; s < neta; s++) for (int t = 0; t < neta; t++) Cee(s, t) = 0.5 * (trace(Hti * d2HtEE[s][t]) - trace(Hti * dHtE[s] * Hti * dHtE[t]));
  // ---- accessors ----
  auto isDir = [&](int p) { return p < ndirP; };
  auto dOf = [&](int p) { return dirP[p] - 1; };
  auto omc = [&](int p) { return p - ndirP; };
  auto McolData = [&](int p) -> vec { if (!isDir(p)) return vec(dOi.slice(omc(p)) * ehat);
    int d = dOf(p); vec r = Nf.col(d);
    for (int l = 0; l < neta; l++) { double v = 0.0; for (int o = 0; o < nobs; o++) v += a(o, l) * (res[o] * iR2[o]) * aRc(o, d); r[l] += v; }
    return r; };
  auto dHtP = [&](int p) -> mat { if (isDir(p)) return dHtDir(dOf(p), aRc); return dOi.slice(omc(p)); };
  auto d2Phi = [&](int aa, int bb) -> double {
    if (!isDir(aa) && !isDir(bb)) return 0.5 * as_scalar(ehat.t() * d2Oi.slice(omc(aa) * nom + omc(bb)) * ehat) + 0.5 * d2LD(omc(aa), omc(bb));
    if (!isDir(aa) || !isDir(bb)) return 0.0;
    int da = dOf(aa), db = dOf(bb); double v = 0.0;
    for (int o = 0; o < nobs; o++) v += rff[o] * a(o, da) * a(o, db) + rf[o] * A(o, da, db) + rfR[o] * (a(o, da) * aRc(o, db) + aRc(o, da) * a(o, db)) +
      rRR[o] * aRc(o, da) * aRc(o, db) + rR[o] * ARc(o, da, db);
    return v; };
  auto SmatEBE = [&](int p) -> mat { mat M(neta, ndir, fill::zeros);
    if (!isDir(p)) { M.cols(0, neta - 1) = dOi.slice(omc(p)); return M; }
    int d = dOf(p);
    for (int l = 0; l < neta; l++) for (int s = 0; s < ndir; s++) { double v = Tnf(l, d, s);
      for (int o = 0; o < nobs; o++) v += -aRc(o, d) * iR2[o] * a(o, s) * a(o, l) + (res[o] * iR2[o]) * aRc(o, d) * A(o, l, s);
      M(l, s) = v; }
    return M; };
  auto SvecEBE = [&](int aa, int bb) -> vec { vec v(neta, fill::zeros);
    bool ta = isDir(aa), tb = isDir(bb);
    if (!ta && !tb) { v = d2Oi.slice(omc(aa) * nom + omc(bb)) * ehat; return v; }
    if (!ta || !tb) return v;
    int da = dOf(aa), db = dOf(bb);
    for (int l = 0; l < neta; l++) { double s = Tnf(l, da, db);
      for (int o = 0; o < nobs; o++) { double w = -(a(o, da) * aRc(o, db) + aRc(o, da) * a(o, db)) * iR2[o] +
          (res[o] * iR2[o]) * ARc(o, da, db) - 2.0 * res[o] * aRc(o, da) * aRc(o, db) * iR3[o];
        s += w * a(o, l) + (res[o] * iR2[o]) * (aRc(o, da) * A(o, l, db) + aRc(o, db) * A(o, l, da)); }
      v[l] = s; }
    return v; };
  auto d2HtEtaP = [&](int p, int l) -> mat { if (!isDir(p)) return zeros<mat>(neta, neta); return d2HtDir(dOf(p), l, aRc, aRe, ARe); };
  auto d2Ht_pp = [&](int aa, int bb) -> mat { bool ta = isDir(aa), tb = isDir(bb);
    if (ta && tb) return d2HtDir(dOf(aa), dOf(bb), aRc, aRc, ARc);
    if (!ta && !tb) return d2Oi.slice(omc(aa) * nom + omc(bb));
    return zeros<mat>(neta, neta); };
  auto Cpe = [&](int p, int l) { return 0.5 * (trace(Hti * d2HtEtaP(p, l)) - trace(Hti * dHtP(p) * Hti * dHtE[l])); };
  auto Cpp = [&](int aa, int bb) { return 0.5 * (trace(Hti * d2Ht_pp(aa, bb)) - trace(Hti * dHtP(aa) * Hti * dHtP(bb))); };
  std::vector<vec> McolV(np); for (int p = 0; p < np; p++) McolV[p] = McolData(p);
  mat etaP(neta, np); for (int p = 0; p < np; p++) etaP.col(p) = -HfInv * McolV[p];
  auto eta2 = [&](int aa, int bb) -> vec {
    mat SmA = SmatEBE(aa).cols(0, neta - 1), SmB = SmatEBE(bb).cols(0, neta - 1);
    vec b = SvecEBE(aa, bb) + SmA * etaP.col(bb) + SmB * etaP.col(aa);
    for (int l = 0; l < neta; l++) { double acc = 0.0;
      for (int s = 0; s < neta; s++) for (int t = 0; t < neta; t++) acc += etaP(s, aa) * Tnf(l, s, t) * etaP(t, bb);
      b[l] += acc; }
    return vec(-HfInv * b); };
  std::vector<vec> CpeRow(np);
  for (int p = 0; p < np; p++) { vec r(neta); for (int l = 0; l < neta; l++) r[l] = Cpe(p, l); CpeRow[p] = r; }
  mat R(np, np, fill::zeros);
  for (int aa = 0; aa < np; aa++) for (int bb = aa; bb < np; bb++) {
    vec e2 = eta2(aa, bb);
    double dat = d2Phi(aa, bb) + dot(McolV[aa], etaP.col(bb)) + dot(McolV[bb], etaP.col(aa)) +
      as_scalar(etaP.col(aa).t() * H * etaP.col(bb)) + dot(gPhi, e2);
    double ld = Cpp(aa, bb) + dot(CpeRow[aa], etaP.col(bb)) + dot(CpeRow[bb], etaP.col(aa)) +
      as_scalar(etaP.col(aa).t() * Cee * etaP.col(bb)) + dot(Cen, e2);
    R(aa, bb) = R(bb, aa) = dat + ld;
  }
  return R;
}
