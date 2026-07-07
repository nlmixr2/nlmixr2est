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
