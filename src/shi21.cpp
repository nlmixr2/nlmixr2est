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
  bool finiteF4 = f4.is_finite();
  f1 = f(tp1, id);
  bool finiteF1 = f1.is_finite();
  if (!finiteF4 && finiteF1) {
    // F4 isn't finite, but F1 is finite
    return -4.0;
  } else if (finiteF4 && !finiteF1) {
    return -1.0;
  } else if (!finiteF4 && !finiteF1) {
    return -14.0;
  }
  arma::vec all = abs(f4-4*f1+3*f0)/(8.0*ef);
  if (all.size() == 1) {
    return all(0);
  }
  // Return harmonic mean
  all = 1.0/all;
  double sum = 0.0;
  for (unsigned int j = all.size(); j--;) {
    sum += all[j];
  }
  return (((double)all.size())/sum);
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
  }
  double h0=h;
  double l = 0, u = R_PosInf, rcur = NA_REAL, tmp;

  arma::vec f1(f0.size());
  double lasth = h;
  int iter=0;
  
  while(true) {
    iter++;
    if (iter > maxiter) {
      break;
    }
    rcur = shiRF(h, f, ef, t, id, idx, f0, f1, l, u);
    lasth = h;
    if (rcur == -4.0) {
      // t + 4*h has problems
      // t + 1*h does not
      // hnew = t + 2.5*hold
      h = 2.5/4*h;
      continue;
    } else if (rcur == -1.0) {
      // t + 4*h does not have problems
      // t + 1*h has problems
      h = 1.5*h;
      continue;
    } else if (rcur == -14.0) {
      // t + 4*h has problems
      // t + 1*h has problems
      h = lasth;
      // Calculate gradient again
      rcur = shiRF(h, f, ef, t, id, idx, f0, f1, l, u);
      break;
    } else if (rcur < rl) {
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
  arma::vec fp3 = f(tp3, id);
  fm1 = f(tm1, id);
  arma::vec fm3 = f(tm3, id);
  finiteFp1 = fp1.is_finite();
  finiteFp3 = fp3.is_finite();
  finiteFm1 = fp1.is_finite();
  finiteFm3 = fp3.is_finite();
  if (finiteFp1 && finiteFp3 &&
      finiteFm3 && finite) {
    arma::vec all = abs(fp3-3*fp1+3*fm1-fm3)/(8.0*ef);
    if (fm3.size() == 1) {
      return all(0);
    }
    // Return harmonic mean
    all = 1.0/all;
    double sum = 0.0;
    for (unsigned int j = all.size(); j--;) {
      sum += all[j];
    }
    return (((double)all.size())/sum);
  } else {
    return -1.0;
 }
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
  if (h == 0) {
    h = pow(ef, 0.3333333333333333333333); 
  }
  double h0=h, tmp = h;
  double l = 0, u = R_PosInf, rcur = NA_REAL;
  double hlast = h;

  arma::vec fp1(f0.size());
  arma::vec fm1(f0.size());

  int iter=0;
  bool finiteFp1 = true, finiteFp3 = true,
    finiteFm1=true, finiteFm3=true;
  while(true) {
    iter++;
    if (iter > maxiter) {
      break;
    }
    rcur = shiRC(h, f, ef, t, id, idx, fp1, fm1, l, u,
                 finiteFp1, finiteFp3, finiteFm1, finiteFm3);
    if (rcur == 1.0) {
      if (!finiteFm1 && !finiteFp1 &&
          !finiteFm3 && !finiteFp3) {
        h = hlast;
        // Calculate gradient again
        rcur = shiRC(h, f, ef, t, id, idx, fp1, fm1, l, u,
                     finiteFp1, finiteFp3, finiteFm1, finiteFm3);
        break;
      }
      if (!finiteFm1 || !finiteFp1) {
        // hnew*3 = hold*0.5
        tmp = h*0.5/3.0;
      } else if (!finiteFm3 || !finiteFp3) {
        // hnew*3 = hold*2
        tmp = h*2.0/3.0;
      }
      continue;
    }  else if (rcur < rl) {
      l = h;
    } else if (rcur > ru) {
      u = h;
    } else {
      break;
    }
    hlast = h;
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
