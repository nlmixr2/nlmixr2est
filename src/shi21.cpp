// https://github.com/cran/Matrix/blob/master/R/nearPD.R
// ADAPTIVE FINITE-DIFFERENCE INTERVAL ESTIMATION FOR NOISY DERIVATIVE-FREE OPTIMIZATION
// HAO-JUN MICHAEL SHI, YUCHEN XIE, MELODY QIMING XUAN AND JORGE NOCEDAL
//
// https://arxiv.org/pdf/2110.06380.pdf

#define ARMA_DONT_PRINT_ERRORS
#define STRICT_R_HEADER
#include "armahead.h"
#include "shi21.h"

double shiRF(double &h, shi21fn_type f, double ef, arma::vec &t, int &id, int &idx,
            arma::vec &f0, arma::vec &f1, double &l, double &u) {

  arma::vec tp4 = t;
  arma::vec tp1 = t;
  tp4(idx) += 4*h;
  tp1(idx) += h;
  arma::vec f4 = f(tp4, id);
  f1 = f(tp1, id);
  return fabs(f4(idx)-4*f1(idx)+3*f0(idx))/(8.0*ef);
}

double shi21Forward(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, int maxiter) {
  // Algorithm 2.1 in paper
  h = nm2divSqrt3*sqrt(ef);
  double h0=h;
  double l = 0, u = R_PosInf, rcur = NA_REAL;

  arma::vec f1(f0.size());
  
  int iter=0;
  while(true) {
    iter++;
    if (iter > maxiter) {
      h = h0;
      break;
    }
    rcur = shiRF(h, f, ef, t, id, idx, f0, f1, l, u);
    if (rcur < rl) {
      l = h;
    } else if (rcur > ru) {
      u = h;
    } else {
      break;
    }
    if (!R_finite(u)) {
      h = 4.0*h;
    } else if (l == 0) {
      h = h/4.0;
    } else {
      h = (l + u)/2.0;
    }
  }
  // Need f1 from shiRF to compute forward difference
  gr = (f1-f0)/h;
  return h;
}

double shiRC(double &h, shi21fn_type f, double ef, arma::vec &t, int &id, int &idx,
             arma::vec &fp1, arma::vec &fm1, double &l, double &u) {
  arma::vec tp3 = t;
  arma::vec tp1 = t;
  arma::vec tm3 = t;
  arma::vec tm1 = t;
  tp3(idx)  += 3*h;
  tp1(idx)  += h;
  tm3(idx)  -= 3*h;
  tm1(idx)  -= h;
  fp1 = f(tp1, id);
  arma::vec fp3 = f(tp3, id);
  fm1 = f(tm1, id);
  arma::vec fm3 = f(tm3, id);
  return fabs(fp3(idx)-3*fp1(idx)+3*fm1(idx)-fm3(idx))/(8.0*ef);  
}

double shi21Central(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, double nu,
                    int maxiter) {
  // Algorithm 3.1
  // weights = -0.5, 0.5
  // s = -1, 1
  h = pow(3*ef, 0.3333333333333333333333); // maybe a different value?
  double h0=h;
  double l = 0, u = R_PosInf, rcur = NA_REAL;

  arma::vec fp1(f0.size());
  arma::vec fm1(f0.size());

  int iter=0;
  while(true) {
    iter++;
    if (iter > maxiter) {
      h = h0;
      break;
    }
    rcur = shiRC(h, f, ef, t, id, idx, fp1, fm1, l, u);
    if (rcur < rl) {
      l = h;
    } else if (rcur > ru) {
      u = h;
    } else {
      break;
    }
    if (!R_finite(u)) {
      h = nu*h;
    } else if (l == 0) {
      h = h/nu;
    } else {
      h = (l + u)/2.0;
    }
  }
  // Need f1 from shiRF to compute forward difference
  gr = (fp1 - fm1)/(2*h);
  return h;
}
