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
#include "censEst.h"   // censNormalPartials: exact censored rho(f,R) partials (M2/M3/M4)
using namespace arma;

// Fill the per-observation determinant coefficients dff/dfr/drr (2nd order) and their
// (f,R) partials p_fff/p_ffR/p_frf/p_fRR/p_RRR (3rd order), plus override the score/realized-H
// partials rf/rR/rff/rfR/rRR for censored observations.  Normal obs keep the Gaussian
// expected-info determinant (dff=1/R, dfr=0, drr=0.5/R^2; p_ffR=-1/R^2, p_RRR=-1/R^3,
// rest 0).  Censored obs always get the censored realized 2nd derivs (rf..rRR); the
// determinant coeffs follow only under censOption "laplace" (1) -- "gauss" (0) keeps the
// Gaussian determinant.  censv[o]!=0 or a finite limv[o] marks a censored (M2/M3/M4) obs.
//
// p_frf = d(dfr)/df is carried separately from p_ffR = d(dff)/dR.  The two agree only when
// (dff,dfr,drr) are second partials of a potential, which the censored laplace determinant is
// (cp[2..4] are exact 2nd derivs of the censored log-density, cp[5..8] their 3rd derivs) and
// the Gauss-Newton expected-info determinant is not: there dfr is identically zero in f and R,
// so d(dfr)/df = 0 while d(dff)/dR = -1/R^2.
static inline void censGradCoefs(const arma::ivec& censv, const arma::vec& limv,
                                 const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                                 int censOpt, int nobs,
                                 arma::vec& rf, arma::vec& rR, arma::vec& rff, arma::vec& rfR, arma::vec& rRR,
                                 arma::vec& dff, arma::vec& dfr, arma::vec& drr,
                                 arma::vec& pfff, arma::vec& pffR, arma::vec& pfRR, arma::vec& pRRR,
                                 arma::vec& pfrf) {
  const bool hasCens = ((int)censv.n_elem == nobs);
  // normal defaults for the determinant coeffs
  dff = 1.0 / Rv;  dfr = arma::zeros<arma::vec>(nobs);  drr = 0.5 / square(Rv);
  pfff = arma::zeros<arma::vec>(nobs);  pffR = -1.0 / square(Rv);
  pfRR = arma::zeros<arma::vec>(nobs);  pRRR = -1.0 / pow(Rv, 3);
  pfrf = arma::zeros<arma::vec>(nobs);   // dfr == 0 in f and R, so d(dfr)/df = 0 (not p_ffR)
  if (!hasCens) return;
  for (int o = 0; o < nobs; o++) {
    double lim = limv.n_elem == (unsigned) nobs ? limv[o] : R_NegInf;
    int cens = censv[o];
    bool isCens = (cens != 0) || (R_FINITE(lim) && !ISNA(lim));
    if (!isCens) continue;
    double cp[9]; for (int i = 0; i < 9; i++) cp[i] = 0.0;
    censNormalPartials((double)cens, yv[o], lim, fv[o], Rv[o], 3, cp);
    rf[o] = cp[0];  rR[o] = cp[1];  rff[o] = cp[2];  rfR[o] = cp[3];  rRR[o] = cp[4];  // always
    if (censOpt == 1) {   // laplace: exact censored determinant
      dff[o] = cp[2];  dfr[o] = cp[3];  drr[o] = cp[4];
      pfff[o] = cp[5]; pffR[o] = cp[6]; pfRR[o] = cp[7]; pRRR[o] = cp[8];
      pfrf[o] = cp[6];   // censored determinant is integrable, so d(dfr)/df = p_ffR
    }
  }
}

// Overwrite the per-obs rho SCORE partials (1st..3rd order) with the exact censored
// (M2/M3/M4) values on censored observations; normal obs keep the Gaussian forms.  Used by
// the (f,R) covariance score terms (Gdd/N/Tn -> the true inner Hessian and its parameter
// chain).  The Laplace-determinant block (Ht/dHtD/d2HtDD) stays Gauss-Newton, matching the
// default censOption="gauss" fit; a laplace-censored cov (censored determinant) bows out to
// FD in R.  censv[o]!=0 or a finite limv[o] marks a censored obs.  rfff is Gaussian-zero.
static inline void censScoreCoefs(const arma::ivec& censv, const arma::vec& limv,
                                  const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv, int nobs,
                                  arma::vec& rf, arma::vec& rR, arma::vec& rff, arma::vec& rfR, arma::vec& rRR,
                                  arma::vec& rffR, arma::vec& rfRR, arma::vec& rRRR, arma::vec& rfff) {
  const bool hasCens = ((int)censv.n_elem == nobs);
  rfff = arma::zeros<arma::vec>(nobs);            // Gaussian rho has rho_fff = 0
  if (!hasCens) return;
  for (int o = 0; o < nobs; o++) {
    double lim = limv.n_elem == (unsigned) nobs ? limv[o] : R_NegInf;
    int cens = censv[o];
    bool isCens = (cens != 0) || (R_FINITE(lim) && !ISNA(lim));
    if (!isCens) continue;
    double cp[9]; for (int i = 0; i < 9; i++) cp[i] = 0.0;
    censNormalPartials((double)cens, yv[o], lim, fv[o], Rv[o], 3, cp);
    rf[o] = cp[0]; rR[o] = cp[1]; rff[o] = cp[2]; rfR[o] = cp[3]; rRR[o] = cp[4];
    rfff[o] = cp[5]; rffR[o] = cp[6]; rfRR[o] = cp[7]; rRRR[o] = cp[8];
  }
}

// FOCE (frozen-R0) censored SCORE: override rho_f, rho_R and the inner-Hessian 2nd derivative
// rff (=1/R0 for a normal obs, exact at (f,R0) for a censored obs) on censored observations.
// The FOCE inner Hessian freezes the variance, so only the f-chain (rho_f, rho_ff) enters the
// eta-block; rho_R feeds the parameter columns.  The determinant stays Gauss-Newton (gauss).
static inline void censFoceScoreCoefs(const arma::ivec& censv, const arma::vec& limv,
                                      const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v, int nobs,
                                      arma::vec& rho_f, arma::vec& rho_R, arma::vec& rff, arma::vec& rfR) {
  rff = 1.0 / R0v;
  rfR = (yv - fv) / arma::square(R0v);                 // normal d(rho_f)/dR0 (the R0-chain cross deriv)
  const bool hasCens = ((int)censv.n_elem == nobs);
  if (!hasCens) return;
  for (int o = 0; o < nobs; o++) {
    double lim = limv.n_elem == (unsigned) nobs ? limv[o] : R_NegInf;
    int cens = censv[o];
    bool isCens = (cens != 0) || (R_FINITE(lim) && !ISNA(lim));
    if (!isCens) continue;
    double cp[9]; for (int i = 0; i < 9; i++) cp[i] = 0.0;
    censNormalPartials((double)cens, yv[o], lim, fv[o], R0v[o], 2, cp);
    rho_f[o] = cp[0]; rho_R[o] = cp[1]; rff[o] = cp[2]; rfR[o] = cp[3];
  }
}

// R-callable wrapper for the exact censored rho(f,R) partials (M2/M3/M4).  Returns an
// nobs x 9 matrix of rho_{f,r,ff,fr,rr,fff,ffr,frr,rrr}; used by the FOCE EBE re-solve to
// build the censored inner score/Hessian (q0=rho_f, q1=rho_ff) at the frozen variance.
// [[Rcpp::export]]
arma::mat censNormalPartials_(const arma::ivec& cens, const arma::vec& dv, const arma::vec& lim,
                              const arma::vec& fv, const arma::vec& rv, int order) {
  const int n = (int)fv.n_elem;
  arma::mat out(n, 9, arma::fill::zeros);
  for (int o = 0; o < n; o++) {
    double l = (lim.n_elem == (unsigned) n) ? lim[o] : R_NegInf;
    double cp[9]; for (int i = 0; i < 9; i++) cp[i] = 0.0;
    censNormalPartials((double)cens[o], dv[o], l, fv[o], rv[o], order, cp);
    for (int i = 0; i < 9; i++) out(o, i) = cp[i];
  }
  return out;
}

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
                                const arma::mat& dvSens,
                                const arma::ivec& censv, const arma::vec& limv, int censOpt,
                                const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                                const arma::vec& ehat, const arma::mat& Oi,
                                const arma::cube& dOiEst, const arma::vec& tr28,
                                int neta, int nth, int nsg, int nom,
                                const arma::ivec& dirTh, const arma::ivec& sigCol,
                                arma::vec& g_out, arma::mat& etaP_out) {
  const int nobs = (int)a.n_rows;
  const int ndir = (int)a.n_cols;
  const int np = nth + nsg + nom;
  // DV-transform chain (estimated boxCox/yeoJohnson lambda): a lambda DIRECTION also
  // moves the DV (y'=tbs(DV,lambda)), so the residual pred sensitivity is a-dvSens
  // (dvSens=dy'/dlambda, nonzero only in the lambda column; the determinant keeps pred a).
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir && dvSens.n_rows == (unsigned) nobs);
  // per-observation rho(f,R,y) partials (normal); censored obs override rf..rRR below
  vec res = yv - fv;
  vec rf = -res / Rv, rR = 0.5 * (1.0 / Rv - square(res) / square(Rv));
  vec rff = 1.0 / Rv, rfR = res / square(Rv), rRR = 0.5 * (-1.0 / square(Rv) + 2.0 * square(res) / pow(Rv, 3));
  // determinant coefficients dff/dfr/drr (+ 3rd-order coeff partials pfff/pffR/pfRR/pRRR);
  // normal = Gauss-Newton expected info, censored+laplace = exact censored 2nd derivative.
  vec dff, dfr, drr, pfff, pffR, pfRR, pRRR, pfrf;
  censGradCoefs(censv, limv, fv, yv, Rv, censOpt, nobs,
                rf, rR, rff, rfR, rRR, dff, dfr, drr, pfff, pffR, pfRR, pRRR, pfrf);
  // exact inner Hessian H (eta x eta), N (eta x dir), determinant Ht
  mat H = Oi, Ht = Oi, N(neta, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) {
    for (int m = 0; m < neta; m++) {
      double sh = 0.0, sht = 0.0;
      for (int o = 0; o < nobs; o++) {
        sh += rff[o] * a(o, l) * a(o, m) + rfR[o] * (a(o, l) * aR(o, m) + aR(o, l) * a(o, m)) +
          rRR[o] * aR(o, l) * aR(o, m) + rf[o] * A(o, l, m) + rR[o] * AR(o, l, m);
        sht += dff[o] * a(o, l) * a(o, m) + dfr[o] * (a(o, l) * aR(o, m) + aR(o, l) * a(o, m)) +
          drr[o] * aR(o, l) * aR(o, m);
      }
      H(l, m) += sh; Ht(l, m) += sht;
    }
    for (int d = 0; d < ndir; d++) {
      double s = 0.0;
      for (int o = 0; o < nobs; o++) {
        double ad = hasDv ? a(o, d) - dvSens(o, d) : a(o, d);  // residual pred sensitivity
        s += rff[o] * a(o, l) * ad + rfR[o] * (a(o, l) * aR(o, d) + aR(o, l) * ad) +
          rRR[o] * aR(o, l) * aR(o, d) + rf[o] * A(o, l, d) + rR[o] * AR(o, l, d);
      }
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
        for (int o = 0; o < nobs; o++) {
          // d(dff)/ddir_s = pfff*a(s) + pffR*aR(s); d(dfr)/ddir_s = pfrf*a(s) + pfRR*aR(s);
          // d(drr)/ddir_s = pfRR*a(s) + pRRR*aR(s)  (coeff (f,R) partials chained through dir s).
          // The f-partial of dfr is pfrf, not pffR -- they differ for the Gauss-Newton
          // determinant, where dfr is identically zero (see censGradCoefs).
          double ddff = pfff[o] * a(o, s) + pffR[o] * aR(o, s);
          double ddfr = pfrf[o] * a(o, s) + pfRR[o] * aR(o, s);
          double ddrr = pfRR[o] * a(o, s) + pRRR[o] * aR(o, s);
          v += ddff * a(o, l) * a(o, m) + dff[o] * (A(o, l, s) * a(o, m) + a(o, l) * A(o, m, s)) +
            ddfr * (a(o, l) * aR(o, m) + aR(o, l) * a(o, m)) +
            dfr[o] * (A(o, l, s) * aR(o, m) + a(o, l) * AR(o, m, s) + AR(o, l, s) * a(o, m) + aR(o, l) * A(o, m, s)) +
            ddrr * aR(o, l) * aR(o, m) + drr[o] * (AR(o, l, s) * aR(o, m) + aR(o, l) * AR(o, m, s));
        }
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
        for (int o = 0; o < nobs; o++) {
          // sigma moves R only: d(dff)/dsig = pffR*Rsig, d(dfr)/dsig = pfRR*Rsig,
          // d(drr)/dsig = pRRR*Rsig; aR moves via RsigDir, a is fixed.
          double ddff = pffR[o] * Rsig(o, c), ddfr = pfRR[o] * Rsig(o, c), ddrr = pRRR[o] * Rsig(o, c);
          v += ddff * a(o, l) * a(o, m) +
            ddfr * (a(o, l) * aR(o, m) + aR(o, l) * a(o, m)) +
            dfr[o] * (a(o, l) * RsigDir(o, m, c) + RsigDir(o, l, c) * a(o, m)) +
            ddrr * aR(o, l) * aR(o, m) + drr[o] * (RsigDir(o, l, c) * aR(o, m) + aR(o, l) * RsigDir(o, m, c));
        }
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
      double s = 0.0; for (int o = 0; o < nobs; o++)
        s += rf[o] * (hasDv ? a(o, d) - dvSens(o, d) : a(o, d)) + rR[o] * aR(o, d);
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
                               const arma::mat& dvSens,
                               const arma::ivec& censv, const arma::vec& limv, int censOpt,
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                               const arma::vec& ehat, const arma::mat& Oi,
                               const arma::cube& dOiEst, const arma::vec& tr28,
                               int neta, int nth, int nsg, int nom,
                               const arma::ivec& dirTh, const arma::ivec& sigCol) {
  vec g; mat etaP;
  foceiGradSubjectFR_(a, A, aR, AR, Rsig, RsigDir, dvSens, censv, limv, censOpt, fv, yv, Rv, ehat, Oi, dOiEst, tr28,
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
                           const arma::mat& dvSens,
                           const arma::ivec& censv, const arma::vec& limv, int censOpt,
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
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir);
  const bool hasCens = ((int)censv.n_elem == (int)fv.n_elem);
  // arma's inv() throws on a singular H/Ht (Makevars.in sets ARMA_DONT_USE_OPENMP but not
  // ARMA_DONT_USE_EXCEPTIONS), and an exception escaping an OpenMP structured block is
  // std::terminate -- the R PROCESS dies.  It is not an R condition, so the driver's
  // tryCatch(..., error=function(e) NULL) cannot intercept it and the FD fallback never
  // runs.  (innerOpt in inner.cpp wraps its region for exactly this reason.)  Catch per
  // subject and poison it: the result carries NaN and the R driver's is.finite() gate
  // degrades to finite differences, which is the intended behaviour.
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < nsub; i++) {
    try {
      int o0 = obsOffset[i], o1 = obsOffset[i + 1] - 1;
      mat ai = a.rows(o0, o1), aRi = aR.rows(o0, o1);
      cube Ai = A.rows(o0, o1), ARi = AR.rows(o0, o1);
      mat Rsigi = hasSig ? mat(Rsig.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      cube RsigDiri = hasSig ? cube(RsigDir.rows(o0, o1)) : cube(o1 - o0 + 1, ndir, 0);
      mat dvi = hasDv ? mat(dvSens.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      ivec censi = hasCens ? ivec(censv.subvec(o0, o1)) : ivec();
      vec limi = hasCens ? vec(limv.subvec(o0, o1)) : vec();
      vec gi; mat etaPi;
      foceiGradSubjectFR_(ai, Ai, aRi, ARi, Rsigi, RsigDiri, dvi, censi, limi, censOpt, fv.subvec(o0, o1), yv.subvec(o0, o1),
                          Rv.subvec(o0, o1), ehat.row(i).t(), Oi, dOiEst, tr28,
                          neta, nth, nsg, nom, dirTh, sigCol, gi, etaPi);
      gmat.col(i) = gi; etaPall.slice(i) = etaPi;
    } catch (...) {
      gmat.col(i).fill(datum::nan); etaPall.slice(i).fill(datum::nan);
    }
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
// Shared core (called from the single-subject export and the batched OpenMP driver).
static void foceiGradSubjectFoceFR_(const arma::mat& a, const arma::cube& A,
                                    const arma::mat& aRe, const arma::mat& aRc,
                                    const arma::mat& R0sig, const arma::mat& dvSens,
                                    const arma::ivec& censv, const arma::vec& limv,
                                    const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                                    const arma::vec& ehat, const arma::mat& Oi,
                                    const arma::cube& dOiEst, const arma::vec& tr28,
                                    int neta, int nth, int nsg, int nom,
                                    const arma::ivec& dirTh, const arma::ivec& sigCol, int fp,
                                    arma::vec& g_out, arma::mat& etaP_out) {
  const int nobs = (int)a.n_rows;
  const int ndir = (int)a.n_cols;
  const int np = nth + nsg + nom;
  // DV-transform chain: a lambda direction also moves y'=tbs(DV,lambda); the residual
  // pred sensitivity is a-dvSens (dvSens=dy'/dlambda, lambda column only).
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir && dvSens.n_rows == (unsigned) nobs);
  vec res = yv - fv;
  vec rho_f = -res / R0v, rho_R = 0.5 * (1.0 / R0v - square(res) / square(R0v));
  // censored (M2/M3/M4) score overrides rho_f/rho_R, the frozen-R inner 2nd deriv rff, and the
  // R0-chain cross deriv rfR = d(rho_f)/dR0 (used by the EBE-sensitivity McolEBE).
  vec rff, rfR;
  censFoceScoreCoefs(censv, limv, fv, yv, R0v, nobs, rho_f, rho_R, rff, rfR);
  vec q0 = rho_f, q1 = rff;
  vec iR = 1.0 / R0v, iR2 = square(iR);
  vec gPhi = Oi * ehat;
  for (int l = 0; l < neta; l++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += rho_f[o] * a(o, l) + rho_R[o] * aRe(o, l); gPhi[l] += s; }
  mat Hf = Oi, Ht = Oi, Nf(neta, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) {
    for (int m = 0; m < neta; m++) { double sh = 0.0, st = 0.0;
      for (int o = 0; o < nobs; o++) { sh += q1[o] * a(o, l) * a(o, m) + q0[o] * A(o, l, m); st += a(o, l) * a(o, m) * iR[o]; }
      Hf(l, m) += sh; Ht(l, m) += st; }
    for (int d = 0; d < ndir; d++) { double s = 0.0; for (int o = 0; o < nobs; o++) {
      double ad = hasDv ? a(o, d) - dvSens(o, d) : a(o, d); s += q1[o] * a(o, l) * ad + q0[o] * A(o, l, d); } Nf(l, d) = s; }
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
  // d(S_FOCE)/dparam explicit R0-chain uses d(rho_f)/dR0 = rfR (censored-exact; = res/R0^2 normal)
  auto McolEBE = [&](int p) -> vec { int t = typ(p);
    if (t == 0) { int d = dirTh[p] - 1; vec r = Nf.col(d);
      for (int l = 0; l < neta; l++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += a(o, l) * rfR[o] * aRc(o, d); r[l] += s; } return r; }
    if (t == 1) { int c = sigCol[p - nth] - 1; vec r(neta);
      for (int l = 0; l < neta; l++) { double s = 0.0; for (int o = 0; o < nobs; o++) s += a(o, l) * rfR[o] * R0sig(o, c); r[l] = s; } return r; }
    return vec(dOiEst.slice(p - nth - nsg) * ehat); };
  auto dHt_p = [&](int p) -> mat { int t = typ(p);
    if (t == 0) { int d = dirTh[p] - 1; mat D = dHtD[d]; if (!fp) D += ouRc(-aRc.col(d) % iR2); return D; }
    if (t == 1) { int c = sigCol[p - nth] - 1; return ouRc(-R0sig.col(c) % iR2); }
    return dOiEst.slice(p - nth - nsg); };
  auto dPhiExplicit = [&](int p) -> double { int t = typ(p);
    if (t == 0) { int d = dirTh[p] - 1; double s = 0.0; for (int o = 0; o < nobs; o++) s += rho_f[o] * (hasDv ? a(o, d) - dvSens(o, d) : a(o, d)) + rho_R[o] * aRc(o, d); return s; }
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
  g_out = g; etaP_out = etaP;
}

// Single-subject export (oracle / R fallback): thin wrapper over foceiGradSubjectFoceFR_.
// [[Rcpp::export]]
Rcpp::List foceiSubjectGradFoceFR_(const arma::mat& a, const arma::cube& A,
                                   const arma::mat& aRe, const arma::mat& aRc,
                                   const arma::mat& R0sig, const arma::mat& dvSens,
                                   const arma::ivec& censv, const arma::vec& limv,
                                   const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                                   const arma::vec& ehat, const arma::mat& Oi,
                                   const arma::cube& dOiEst, const arma::vec& tr28,
                                   int neta, int nth, int nsg, int nom,
                                   const arma::ivec& dirTh, const arma::ivec& sigCol, int fp) {
  vec g; mat etaP;
  foceiGradSubjectFoceFR_(a, A, aRe, aRc, R0sig, dvSens, censv, limv, fv, yv, R0v, ehat, Oi, dOiEst, tr28,
                          neta, nth, nsg, nom, dirTh, sigCol, fp, g, etaP);
  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaP);
}

// Batched (f,R) FOCE outer gradient over ALL subjects in one OpenMP-parallel C++ call.
// aRe/aRc/R0sig are the per-subject frozen-R0 sensitivities resolved in R (from E/E0),
// concatenated over observations (obsOffset[i]..obsOffset[i+1]-1 are subject i's rows);
// ehat is nsub x neta.  Returns the summed gradient g (np) and per-subject etaP cube.
// [[Rcpp::export]]
Rcpp::List foceiGradAllFoceFR_(const arma::mat& a, const arma::cube& A,
                               const arma::mat& aRe, const arma::mat& aRc, const arma::mat& R0sig,
                               const arma::mat& dvSens, const arma::ivec& censv, const arma::vec& limv,
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                               const arma::mat& ehat, const arma::ivec& obsOffset,
                               const arma::mat& Oi, const arma::cube& dOiEst, const arma::vec& tr28,
                               int neta, int nth, int nsg, int nom,
                               const arma::ivec& dirTh, const arma::ivec& sigCol, int fp, int ncores) {
  const int nsub = (int)ehat.n_rows;
  const int np = nth + nsg + nom;
  const int ndir = (int)a.n_cols;
  mat gmat(np, nsub, fill::zeros);
  cube etaPall(neta, np, nsub, fill::zeros);
  const bool hasSig = (R0sig.n_cols > 0);
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir);
  const bool hasCens = ((int)censv.n_elem == (int)fv.n_elem);
  // arma's inv() throws on a singular H/Ht (Makevars.in sets ARMA_DONT_USE_OPENMP but not
  // ARMA_DONT_USE_EXCEPTIONS), and an exception escaping an OpenMP structured block is
  // std::terminate -- the R PROCESS dies.  It is not an R condition, so the driver's
  // tryCatch(..., error=function(e) NULL) cannot intercept it and the FD fallback never
  // runs.  (innerOpt in inner.cpp wraps its region for exactly this reason.)  Catch per
  // subject and poison it: the result carries NaN and the R driver's is.finite() gate
  // degrades to finite differences, which is the intended behaviour.
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < nsub; i++) {
    try {
      int o0 = obsOffset[i], o1 = obsOffset[i + 1] - 1;
      mat ai = a.rows(o0, o1), aRei = aRe.rows(o0, o1), aRci = aRc.rows(o0, o1);
      cube Ai = A.rows(o0, o1);
      mat R0sigi = hasSig ? mat(R0sig.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      mat dvi = hasDv ? mat(dvSens.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      ivec censi = hasCens ? ivec(censv.subvec(o0, o1)) : ivec();
      vec limi = hasCens ? vec(limv.subvec(o0, o1)) : vec();
      vec gi; mat etaPi;
      foceiGradSubjectFoceFR_(ai, Ai, aRei, aRci, R0sigi, dvi, censi, limi, fv.subvec(o0, o1), yv.subvec(o0, o1),
                              R0v.subvec(o0, o1), ehat.row(i).t(), Oi, dOiEst, tr28,
                              neta, nth, nsg, nom, dirTh, sigCol, fp, gi, etaPi);
      gmat.col(i) = gi; etaPall.slice(i) = etaPi;
    } catch (...) {
      gmat.col(i).fill(datum::nan); etaPall.slice(i).fill(datum::nan);
    }
  }
  vec g = sum(gmat, 1);
  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaPall);
}

// (f,R) FOCEI per-subject observed-information R (oracle: .foceiAnalyticSubjectRFR).
// Analytic 1st/2nd-order sensitivities a/A (prediction) + aR/AR (variance); the 3rd-order
// tensors Ath/AthR come from Shi-FD and are passed reshaped to cubes (nobs, ndir, ndir*ndir):
// Ath[o,l,s,t] == Ath(o, l, s + t*ndir).  Every non-Omega param is a direction (dirP,
// 1-based); a residual sigma direction has a=A=Ath=0.  Omega derivatives: dOi
// (neta,neta,nom), d2Oi (neta,neta,nom*nom) with slice a*nom+b, d2LD (nom,nom).
// Shared core (called from the single-subject export and the batched OpenMP driver).
static arma::mat foceiRSubjectFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                           const arma::mat& aR, const arma::cube& AR, const arma::cube& AthR,
                           const arma::mat& dvSens, const arma::mat& dvSens2,
                           const arma::ivec& censv, const arma::vec& limv,
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
  // censored (M2/M3/M4) SCORE partials overwrite the Gaussian rho derivs (rfff is 0 for a
  // normal obs); the determinant below stays Gauss-Newton (censOption="gauss").
  vec rfff;
  censScoreCoefs(censv, limv, fv, yv, Rv, nobs, rf, rR, rff, rfR, rRR, rffR, rfRR, rRRR, rfff);
  vec iR = 1.0 / Rv, iR2 = square(iR), iR3 = pow(iR, 3), iR4 = pow(iR, 4);
  auto Ai = [&](const arma::cube& T, int o, int l, int s, int t) { return T(o, l, s + t * ndir); };
  // DV-transform chain (estimated boxCox/yeoJohnson lambda): the residual pred sensitivity
  // ra = a - dy'/dlambda (dvSens, lambda column only) enters the RHO residual terms; the
  // lambda-lambda 2nd derivative also carries d2y'/dlambda2 (dvSens2).  The determinant
  // (Gauss-Newton, no residual) and the eta-block keep the pure pred a.  All-zero (no-op)
  // for non-transform models; the 3rd-order Ath/AthR are DV-free (the DV has no eta chain).
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir && dvSens.n_rows == (unsigned) nobs);
  auto ra = [&](int o, int d) { return hasDv ? a(o, d) - dvSens(o, d) : a(o, d); };
  auto dvY = [&](int o, int d) { return hasDv ? dvSens2(o, d) : 0.0; };  // d2y'/dlambda2 (lambda col)
  // Gdd: (f,R) 2nd total derivative of the density between two directions
  auto Gdd = [&](int da, int db) {
    double s = 0.0;
    for (int o = 0; o < nobs; o++) {
      s += rff[o] * ra(o, da) * ra(o, db) + rfR[o] * (ra(o, da) * aR(o, db) + aR(o, da) * ra(o, db)) +
        rRR[o] * aR(o, da) * aR(o, db) + rf[o] * A(o, da, db) + rR[o] * AR(o, da, db);
      if (da == db) s += -rf[o] * dvY(o, da);         // (r/R) d2y'/dlambda2, lambda-lambda only
    }
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
      double ras = ra(o, s), rat = ra(o, t), Yst = (s == t) ? dvY(o, s) : 0.0;  // DV residual sens + d2y'/dl2
      double us = rff[o] * ras + rfR[o] * aR(o, s), ut = rff[o] * rat + rfR[o] * aR(o, t);
      double ust = rfff[o] * ras * rat + rffR[o] * (ras * aR(o, t) + aR(o, s) * rat) + rfRR[o] * aR(o, s) * aR(o, t) +
        rff[o] * A(o, s, t) - rff[o] * Yst + rfR[o] * AR(o, s, t);
      double ws = rfR[o] * ras + rRR[o] * aR(o, s), wt = rfR[o] * rat + rRR[o] * aR(o, t);
      double wst = rffR[o] * ras * rat + rfRR[o] * (ras * aR(o, t) + aR(o, s) * rat) +
        rRRR[o] * aR(o, s) * aR(o, t) + rfR[o] * A(o, s, t) - rfR[o] * Yst + rRR[o] * AR(o, s, t);
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

// Single-subject export (oracle / R fallback): thin wrapper over foceiRSubjectFR_.
// [[Rcpp::export]]
arma::mat foceiSubjectRFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                           const arma::mat& aR, const arma::cube& AR, const arma::cube& AthR,
                           const arma::mat& dvSens, const arma::mat& dvSens2,
                           const arma::ivec& censv, const arma::vec& limv,
                           const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                           const arma::vec& ehat, const arma::mat& Oi,
                           const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                           int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP) {
  return foceiRSubjectFR_(a, A, Ath, aR, AR, AthR, dvSens, dvSens2, censv, limv, fv, yv, Rv, ehat, Oi, dOi, d2Oi, d2LD,
                          neta, ndir, ndirP, nom, dirP);
}

// Batched (f,R) FOCEI observed-information R summed over ALL subjects in one OpenMP call.
// Sensitivities (incl. the 3rd-order Ath/AthR, already Shi-FD'd per subject) are concatenated
// over observations (obsOffset[i]..obsOffset[i+1]-1 are subject i's rows); ehat is nsub x neta.
// [[Rcpp::export]]
arma::mat foceiRAllFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                       const arma::mat& aR, const arma::cube& AR, const arma::cube& AthR,
                       const arma::mat& dvSens, const arma::mat& dvSens2,
                       const arma::ivec& censv, const arma::vec& limv,
                       const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv,
                       const arma::mat& ehat, const arma::ivec& obsOffset,
                       const arma::mat& Oi, const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                       int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP, int ncores) {
  const int nsub = (int)ehat.n_rows;
  const int np = ndirP + nom;
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir);
  const bool hasCens = ((int)censv.n_elem == (int)fv.n_elem);
  cube Rall(np, np, nsub, fill::zeros);
  // arma's inv() throws on a singular H/Ht (Makevars.in sets ARMA_DONT_USE_OPENMP but not
  // ARMA_DONT_USE_EXCEPTIONS), and an exception escaping an OpenMP structured block is
  // std::terminate -- the R PROCESS dies.  It is not an R condition, so the driver's
  // tryCatch(..., error=function(e) NULL) cannot intercept it and the FD fallback never
  // runs.  (innerOpt in inner.cpp wraps its region for exactly this reason.)  Catch per
  // subject and poison it: the result carries NaN and the R driver's is.finite() gate
  // degrades to finite differences, which is the intended behaviour.
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < nsub; i++) {
    try {
      int o0 = obsOffset[i], o1 = obsOffset[i + 1] - 1;
      mat dvi = hasDv ? mat(dvSens.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      mat dv2i = hasDv ? mat(dvSens2.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      ivec censi = hasCens ? ivec(censv.subvec(o0, o1)) : ivec();
      vec limi = hasCens ? vec(limv.subvec(o0, o1)) : vec();
      Rall.slice(i) = foceiRSubjectFR_(a.rows(o0, o1), A.rows(o0, o1), Ath.rows(o0, o1),
                                       aR.rows(o0, o1), AR.rows(o0, o1), AthR.rows(o0, o1), dvi, dv2i, censi, limi,
                                       fv.subvec(o0, o1), yv.subvec(o0, o1), Rv.subvec(o0, o1),
                                       ehat.row(i).t(), Oi, dOi, d2Oi, d2LD,
                                       neta, ndir, ndirP, nom, dirP);
    } catch (...) {
      Rall.slice(i).fill(datum::nan);
    }
  }
  mat R = sum(Rall, 2);
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
// Shared core (called from the single-subject export and the batched OpenMP driver).
static arma::mat foceiRSubjectFoceFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                               const arma::mat& aRe, const arma::mat& aRc,
                               const arma::cube& ARe, const arma::cube& ARc, const arma::mat& dvSens,
                               const arma::mat& dvSens2,
                               const arma::ivec& censv, const arma::vec& limv,
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                               const arma::vec& ehat, const arma::mat& Oi,
                               const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                               int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP) {
  const int nobs = (int)a.n_rows;
  const int np = ndirP + nom;
  vec res = yv - fv;
  vec rf = -res / R0v, rR = 0.5 * (1.0 / R0v - square(res) / square(R0v));
  vec rff = 1.0 / R0v, rfR = res / square(R0v), rRR = 0.5 * (-1.0 / square(R0v) + 2.0 * square(res) / pow(R0v, 3));
  // R0-chain partials (score side): rffR = d(rho_ff)/dR0, rfRR = d(rho_fR)/dR0.  Normal forms
  // below; censored (M2/M3/M4) obs get the exact values (censScoreCoefs).  rRRR is unused by the
  // frozen-R FOCE cov; rfff feeds the Tnf 3rd-order term (0 for a normal frozen-R obs).
  vec rffR = -1.0 / square(R0v), rfRR = -2.0 * res / pow(R0v, 3), rRRR = arma::zeros<arma::vec>(nobs), rfff;
  censScoreCoefs(censv, limv, fv, yv, R0v, nobs, rf, rR, rff, rfR, rRR, rffR, rfRR, rRRR, rfff);
  vec q0 = rf, q1 = rff;                              // interaction-free inner (censored-aware)
  vec iR = 1.0 / R0v, iR2 = square(iR), iR3 = pow(iR, 3);
  auto Ai = [&](const arma::cube& T, int o, int l, int s, int t) { return T(o, l, s + t * ndir); };
  // DV-transform chain (estimated boxCox/yeoJohnson lambda): the residual pred sensitivity
  // ra = a - dy'/dlambda enters the RHO/inner data terms (Nf, Tnf, d2Phi, Smat/SvecEBE); the
  // lambda-lambda 2nd derivative carries d2y'/dlambda2.  The frozen-R0 determinant and the
  // 3rd-order Ath stay pure.  All-zero (no-op) for non-transform models.  Same DV correction
  // for nonmem-FOCE and foce+ (the DV moves the residual identically; only R0/aRe/aRc differ).
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir && dvSens.n_rows == (unsigned) nobs);
  auto ra = [&](int o, int d) { return hasDv ? a(o, d) - dvSens(o, d) : a(o, d); };
  auto dvY = [&](int o, int d) { return hasDv ? dvSens2(o, d) : 0.0; };  // d2y'/dlambda2 (lambda col)
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
    for (int d = 0; d < ndir; d++) { double v = 0.0; for (int o = 0; o < nobs; o++) v += q1[o] * a(o, l) * ra(o, d) + q0[o] * A(o, l, d); Nf(l, d) = v; } }
  mat HfInv = inv(Hf);
  cube Tnf(neta, ndir, ndir, fill::zeros);
  for (int l = 0; l < neta; l++) for (int s = 0; s < ndir; s++) for (int t = 0; t < ndir; t++) { double v = 0.0;
    for (int o = 0; o < nobs; o++) { double Yst = (s == t) ? dvY(o, s) : 0.0;
      v += rfff[o] * ra(o, s) * ra(o, t) * a(o, l) +   // frozen-R 3rd-order f-chain (0 for normal)
        q1[o] * (A(o, l, s) * ra(o, t) + A(o, l, t) * ra(o, s) + A(o, s, t) * a(o, l) - Yst * a(o, l)) + q0[o] * Ai(Ath, o, l, s, t); }
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
    for (int l = 0; l < neta; l++) { double v = 0.0; for (int o = 0; o < nobs; o++) v += a(o, l) * rfR[o] * aRc(o, d); r[l] += v; }
    return r; };
  auto dHtP = [&](int p) -> mat { if (isDir(p)) return dHtDir(dOf(p), aRc); return dOi.slice(omc(p)); };
  auto d2Phi = [&](int aa, int bb) -> double {
    if (!isDir(aa) && !isDir(bb)) return 0.5 * as_scalar(ehat.t() * d2Oi.slice(omc(aa) * nom + omc(bb)) * ehat) + 0.5 * d2LD(omc(aa), omc(bb));
    if (!isDir(aa) || !isDir(bb)) return 0.0;
    int da = dOf(aa), db = dOf(bb); double v = 0.0;
    for (int o = 0; o < nobs; o++) { double rada = ra(o, da), radb = ra(o, db);
      v += rff[o] * rada * radb + rf[o] * A(o, da, db) + rfR[o] * (rada * aRc(o, db) + aRc(o, da) * radb) +
        rRR[o] * aRc(o, da) * aRc(o, db) + rR[o] * ARc(o, da, db);
      if (da == db) v += -rf[o] * dvY(o, da); }        // (r/R0) d2y'/dlambda2, lambda-lambda only
    return v; };
  auto SmatEBE = [&](int p) -> mat { mat M(neta, ndir, fill::zeros);
    if (!isDir(p)) { M.cols(0, neta - 1) = dOi.slice(omc(p)); return M; }
    int d = dOf(p);
    for (int l = 0; l < neta; l++) for (int s = 0; s < ndir; s++) { double v = Tnf(l, d, s);
      for (int o = 0; o < nobs; o++) v += rffR[o] * aRc(o, d) * ra(o, s) * a(o, l) + rfR[o] * aRc(o, d) * A(o, l, s);
      M(l, s) = v; }
    return M; };
  auto SvecEBE = [&](int aa, int bb) -> vec { vec v(neta, fill::zeros);
    bool ta = isDir(aa), tb = isDir(bb);
    if (!ta && !tb) { v = d2Oi.slice(omc(aa) * nom + omc(bb)) * ehat; return v; }
    if (!ta || !tb) return v;
    int da = dOf(aa), db = dOf(bb);
    for (int l = 0; l < neta; l++) { double s = Tnf(l, da, db);
      for (int o = 0; o < nobs; o++) { double w = rffR[o] * (ra(o, da) * aRc(o, db) + aRc(o, da) * ra(o, db)) +
          rfR[o] * ARc(o, da, db) + rfRR[o] * aRc(o, da) * aRc(o, db);
        s += w * a(o, l) + rfR[o] * (aRc(o, da) * A(o, l, db) + aRc(o, db) * A(o, l, da)); }
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

// Single-subject export (oracle / R fallback): thin wrapper over foceiRSubjectFoceFR_.
// [[Rcpp::export]]
arma::mat foceiSubjectRfoceFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                               const arma::mat& aRe, const arma::mat& aRc,
                               const arma::cube& ARe, const arma::cube& ARc,
                               const arma::mat& dvSens, const arma::mat& dvSens2,
                               const arma::ivec& censv, const arma::vec& limv,
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                               const arma::vec& ehat, const arma::mat& Oi,
                               const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                               int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP) {
  return foceiRSubjectFoceFR_(a, A, Ath, aRe, aRc, ARe, ARc, dvSens, dvSens2, censv, limv, fv, yv, R0v, ehat, Oi, dOi, d2Oi, d2LD,
                              neta, ndir, ndirP, nom, dirP);
}

// Batched (f,R) FOCE observed-information R summed over ALL subjects in one OpenMP call.
// aRe/aRc/ARe/ARc are the per-subject frozen-R0 sensitivities resolved in R (from E/E0),
// concatenated over observations (obsOffset[i]..obsOffset[i+1]-1 are subject i's rows).
// [[Rcpp::export]]
arma::mat foceiRAllFoceFR_(const arma::mat& a, const arma::cube& A, const arma::cube& Ath,
                           const arma::mat& aRe, const arma::mat& aRc,
                           const arma::cube& ARe, const arma::cube& ARc, const arma::mat& dvSens,
                           const arma::mat& dvSens2, const arma::ivec& censv, const arma::vec& limv,
                           const arma::vec& fv, const arma::vec& yv, const arma::vec& R0v,
                           const arma::mat& ehat, const arma::ivec& obsOffset,
                           const arma::mat& Oi, const arma::cube& dOi, const arma::cube& d2Oi, const arma::mat& d2LD,
                           int neta, int ndir, int ndirP, int nom, const arma::ivec& dirP, int ncores) {
  const int nsub = (int)ehat.n_rows;
  const int np = ndirP + nom;
  const bool hasDv = (dvSens.n_cols == (unsigned) ndir);
  const bool hasCens = ((int)censv.n_elem == (int)fv.n_elem);
  cube Rall(np, np, nsub, fill::zeros);
  // arma's inv() throws on a singular H/Ht (Makevars.in sets ARMA_DONT_USE_OPENMP but not
  // ARMA_DONT_USE_EXCEPTIONS), and an exception escaping an OpenMP structured block is
  // std::terminate -- the R PROCESS dies.  It is not an R condition, so the driver's
  // tryCatch(..., error=function(e) NULL) cannot intercept it and the FD fallback never
  // runs.  (innerOpt in inner.cpp wraps its region for exactly this reason.)  Catch per
  // subject and poison it: the result carries NaN and the R driver's is.finite() gate
  // degrades to finite differences, which is the intended behaviour.
#pragma omp parallel for num_threads(ncores)
  for (int i = 0; i < nsub; i++) {
    try {
      int o0 = obsOffset[i], o1 = obsOffset[i + 1] - 1;
      mat dvi = hasDv ? mat(dvSens.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      mat dv2i = hasDv ? mat(dvSens2.rows(o0, o1)) : mat(o1 - o0 + 1, 0);
      ivec censi = hasCens ? ivec(censv.subvec(o0, o1)) : ivec();
      vec limi = hasCens ? vec(limv.subvec(o0, o1)) : vec();
      Rall.slice(i) = foceiRSubjectFoceFR_(a.rows(o0, o1), A.rows(o0, o1), Ath.rows(o0, o1),
                                           aRe.rows(o0, o1), aRc.rows(o0, o1), ARe.rows(o0, o1), ARc.rows(o0, o1), dvi, dv2i,
                                           censi, limi, fv.subvec(o0, o1), yv.subvec(o0, o1), R0v.subvec(o0, o1),
                                           ehat.row(i).t(), Oi, dOi, d2Oi, d2LD,
                                           neta, ndir, ndirP, nom, dirP);
    } catch (...) {
      Rall.slice(i).fill(datum::nan);
    }
  }
  mat R = sum(Rall, 2);
  return R;
}
