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
// [[Rcpp::export]]
Rcpp::List foceiSubjectGradFR_(const arma::mat& a,        // nobs x ndir (d f / d dir)
                               const arma::cube& A,        // nobs x ndir x ndir
                               const arma::mat& aR,        // nobs x ndir (d R / d dir)
                               const arma::cube& AR,       // nobs x ndir x ndir
                               const arma::mat& Rsig,      // nobs x nsg (d R / d sigma)
                               const arma::cube& RsigDir,  // nobs x ndir x nsg (d2 R / d dir d sigma)
                               const arma::vec& fv, const arma::vec& yv, const arma::vec& Rv, // nobs
                               const arma::vec& ehat,      // neta
                               const arma::mat& Oi,        // neta x neta
                               const arma::cube& dOiEst,   // neta x neta x nom
                               const arma::vec& tr28,      // nom
                               int neta, int nth, int nsg, int nom,
                               const arma::ivec& dirTh,    // nth (1-based direction per theta)
                               const arma::ivec& sigCol) { // nsg (1-based Rsig column per sigma)
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
  return Rcpp::List::create(Rcpp::Named("g") = g, Rcpp::Named("etaP") = etaP);
}
