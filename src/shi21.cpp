// ADAPTIVE FINITE-DIFFERENCE INTERVAL ESTIMATION FOR NOISY DERIVATIVE-FREE OPTIMIZATION
// HAO-JUN MICHAEL SHI, YUCHEN XIE, MELODY QIMING XUAN AND JORGE NOCEDAL
//
// https://arxiv.org/pdf/2110.06380.pdf

// Harmonic mean correction see:
// https://github.com/cran/lmomco/blob/8558903cdcdf6ce0640822a8f6ee7caf07ebd451/R/harmonic.mean.R
// https://rdrr.io/cran/lmomco/man/harmonic.mean.html
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "shi21.h"

double shiRF(double &h, shi21fn_type f, double ef, arma::vec &t, int &id, int &idx,
             arma::vec &f0, arma::vec &f1, double &l, double &u,
             bool &finiteF1, bool &finiteF4) {
  arma::vec tp4 = t;
  arma::vec tp1 = t;
  tp4(idx) += 4*h;
  tp1(idx) += h;
  f1 = f(tp1, id);
  finiteF1 = f1.is_finite();
  if (!finiteF1) {
    finiteF4 = true;
    return -1.0;
  }
  arma::vec f4 = f(tp4, id);
  finiteF4 = f4.is_finite();
  if (!finiteF4) {
    return -1.0;
  }
  arma::vec all = abs(f4-4*f1+3*f0)/(8.0*ef);
  if (all.size() == 1) {
    return all(0);
  }
  // Return harmonic mean
  arma::vec all0 = all;
  all = 1.0/all;
  double sum = 0.0;
  int nzero = 0;
  int n = 0;
  for (unsigned int j = all.size(); j--;) {
    if  (all0[j] == 0) {
      nzero++;
    } else {
      sum += all[j];
      n++;
    }
  }
  double correction = (double)(n-nzero)/((double)n);
  if (correction <= 0) correction=1;
  double hm = (double)(n)/sum * correction;
  return hm;
}

double shi21Forward(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, int maxiter) {
  // Algorithm 2.1 in paper
  // q=2, alpha=4, r=3
  // s = 0, 1
  // w = -1, 1
  if (h == 0) {
    h = nm2divSqrt3*sqrt(ef);
  } else {
    h = fabs(h);
  }
  double h0=h;
  double l = 0, u = R_PosInf, rcur = NA_REAL, tmp;
  arma::vec f1(f0.size());
  double lasth = h;
  int iter=0;
  bool finiteF1 = true, finiteF4 = true, calcGrad = false;
  while(true) {
    iter++;
    if (iter > maxiter) {
      h = lasth;
      break;
    }
    rcur = shiRF(h, f, ef, t, id, idx, f0, f1, l, u,
                 finiteF1, finiteF4);
    if (rcur == -1) {
      if (!finiteF1) {
        // hnew = t + 2.5*hold
        h = 0.5*h;
        continue;
      }
      h = 3.5*h;
      if (!calcGrad) {
        lasth = h;
        gr = (f1-f0)/h;
      }
      continue;
    } else {
      lasth = h;
      gr = (f1-f0)/h;
    }
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
  return h;
}

double shiRC(double &h, shi21fn_type f, double ef, arma::vec &t, int &id, int &idx,
             arma::vec &fp1, arma::vec &fm1, double &l, double &u,
             bool &finiteFp1, bool &finiteFp3,
             bool &finiteFm1, bool &finiteFm3) {
  arma::vec tp3 = t;
  arma::vec tp1 = t;
  arma::vec tm3 = t;
  arma::vec tm1 = t;
  tp3(idx)  += 3*h;
  tp1(idx)  += h;
  tm3(idx)  -= 3*h;
  tm1(idx)  -= h;
  fp1 = f(tp1, id);
  finiteFp1 = fp1.is_finite();
  if (!finiteFp1) {
    finiteFm1 = true;
    finiteFp3 = true;
    finiteFm3 = true;
    return -1.0;
  }
  fm1 = f(tm1, id);
  finiteFm1 = fp1.is_finite();
  if (!finiteFm1) {
    finiteFp3 = true;
    finiteFm3 = true;
    return -1.0;
  }
  arma::vec fp3 = f(tp3, id);
  finiteFp3 = fp3.is_finite();
  if (!finiteFp3) {
    finiteFp3 = true;
    return -1.0;
  }
  arma::vec fm3 = f(tm3, id);
  finiteFm3 = fp3.is_finite();
  if (!finiteFm3) {
    return -1.0;
  }
  arma::vec all = abs(fp3-3*fp1+3*fm1-fm3)/(8.0*ef);
  if (fm3.size() == 1) {
    return all(0);
  }
  // Return harmonic mean
  arma::vec all0 = all;
  all = 1.0/all;
  double sum = 0.0;
  int nzero = 0;
  int n = 0;
  for (unsigned int j = all.size(); j--;) {
    if  (all0[j] == 0) {
      nzero++;
    } else {
      sum += all[j];
      n++;
    }
  }
  double correction = (double)(n-nzero)/((double)n);
  if (correction <= 0) correction=1;
  double hm = (double)(n)/sum * correction;
  return hm;
}

double shi21Central(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, double nu,
                    int maxiter) {
  // Algorithm 3.1
  // weights = -0.5, 0.5
  // s = -1, 1
  // Equation 3.3
  //
  if (h == 0.0) {
    h = pow(3.0*ef, 0.3333333333333333333333); 
  } else {
    h = fabs(h);
  }
  double h0=h, tmp = h;
  double l = 0, u = R_PosInf, rcur = NA_REAL;
  double hlast = h;

  arma::vec fp1(f0.size());
  arma::vec fm1(f0.size());

  int iter=0;
  bool finiteFp1 = true, finiteFp3 = true,
    finiteFm1=true, finiteFm3=true, calcGrad=false;
  while(true) {
    iter++;
    if (iter > maxiter) {
      h=hlast;
      break;
    }
    rcur = shiRC(h, f, ef, t, id, idx, fp1, fm1, l, u,
                 finiteFp1, finiteFp3, finiteFm1, finiteFm3);
    // Need f1 from shiRF to compute forward difference
    if (rcur == -1.0) {
      if (!finiteFp1) {
        // hnew*3 = hold*0.5
        h = h*0.5/3.0;
        continue;
      } else if (!finiteFm1) {
        if (!calcGrad) {
          // forward difference
          calcGrad = true;
          gr = (fp1-f0)/h;
        }
        h = h*0.5/3.0;
        continue;
      }
      // hnew*3 = hold*2
      h = h*2.0/3.0;
      if (!calcGrad) {
        // central difference
        calcGrad = true;
        gr = (fp1-fm1)/(2*h);
        hlast = h;
      }
      continue;
    } else {
      calcGrad = true;
      gr = (fp1-fm1)/(2*h);
      hlast = h;      
    }
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
  return h;
}
