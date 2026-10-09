// ADAPTIVE FINITE-DIFFERENCE INTERVAL ESTIMATION FOR NOISY DERIVATIVE-FREE OPTIMIZATION
// HAO-JUN MICHAEL SHI, YUCHEN XIE, MELODY QIMING XUAN AND JORGE NOCEDAL
//
// https://arxiv.org/pdf/2110.06380.pdf

// Components with r < 1 are non-detects: the third difference is within the 8*ef
// noise bound, so r is left-censored at 1 (#1188).  shiRatioCensor picks how they
// enter the harmonic mean:
//   0 "current"    harmonic mean of the nonzero ratios, legacy zero correction
//   1 "detected"   harmonic mean of the detected (r >= 1) ratios; max(r) if none
//   2 "substitute" censored ratios set to the detection limit 1
//   3 "lmomco"     detected harmonic mean times (N - N0)/N, N0 = censored count
// lmomco::harmonic.mean() (TCEQ RG-194 / EPA DFLOW) is the source of option 3.
#define ARMA_WARN_LEVEL 1
#define STRICT_R_HEADER
#include "armahead.h"
#include "shi21.h"

static int shiRatioCensor_ = 0;

// Selects the shiRatio() treatment of censored ratios; returns the previous one.
//[[Rcpp::export]]
int shi21RatioCensorSet(int type) {
  if (type < 0 || type > 3) {
    Rcpp::stop("unknown shi21 ratio censor type");
  }
  int old = shiRatioCensor_;
  shiRatioCensor_ = type;
  return old;
}

static double shiRatioType(const arma::vec &all, int type) {
  if (all.size() == 1) {
    return all(0);
  }
  double sum = 0.0;
  int nzero = 0;
  int n = 0;
  if (type == 0) {
    for (unsigned int j = all.size(); j--;) {
      if  (all[j] == 0) {
        nzero++;
      } else {
        sum += 1.0/all[j];
        n++;
      }
    }
    double correction = (double)(n-nzero)/((double)n);
    if (correction <= 0) correction=1;
    return (double)(n)/sum * correction;
  }
  double rmax = 0.0;
  for (unsigned int j = all.size(); j--;) {
    if (all[j] >= 1.0) {
      sum += 1.0/all[j];
      n++;
    } else {
      nzero++;
      if (all[j] > rmax) rmax = all[j];
      if (type == 2) sum += 1.0;
    }
  }
  switch (type) {
  case 1:
    if (n == 0) return rmax;
    return (double)(n)/sum;
  case 2:
    return (double)(all.size())/sum;
  default:
    if (n == 0) return 0.0;
    return (double)(n)/sum * (double)(n)/((double)(all.size()));
  }
}

// The ratio test statistic of a step: the one ratio itself for a scalar function,
// else a harmonic mean of the per-element ratios (see shiRatioCensor_).
static double shiRatio(const arma::vec &all) {
  return shiRatioType(all, shiRatioCensor_);
}

// R-callable shiRatio() for unit tests.
//[[Rcpp::export]]
double shi21RatioTest(arma::vec all, int type) {
  return shiRatioType(all, type);
}

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
  return shiRatio(abs(f4-4*f1+3*f0)/(8.0*ef));
}

double shi21Forward(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, int maxiter,
                    double hMax, double hMin) {
  // Algorithm 2.1 in paper
  // q=2, alpha=4, r=3
  // s = 0, 1
  // w = -1, 1
  if (h == 0) {
    h = nm2divSqrt3*sqrt(ef);
  } else {
    h = fabs(h);
  }
  // Bound the FD step both ways (hMax/hMin, caller-configurable) -- see shi21Central;
  // an unbounded step corrupts the shared solver state via a degenerate probe, and a
  // vanishing step is roundoff.
  double l = 0, u = R_PosInf, rcur = NA_REAL;
  arma::vec f1(f0.size());
  gr.zeros(f0.n_elem); // avoid uninitialized read if no finite forward step is found
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
      if (h >= hMax) break;
      h = 4.0*h;
      if (h > hMax) h = hMax;
    } else if (l == 0) {
      if (h <= hMin) break;
      h = h/4.0;
      if (h < hMin) h = hMin;
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
  finiteFm1 = fm1.is_finite();
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
  finiteFm3 = fm3.is_finite();
  if (!finiteFm3) {
    return -1.0;
  }
  return shiRatio(abs(fp3-3*fp1+3*fm1-fm3)/(8.0*ef));
}

double shi21Central(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, double nu,
                    int maxiter, double hMax, double hMin) {
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
  // Bound the finite-difference step to a reasonable region (hMax/hMin, caller-
  // configurable).  Too big: a flat objective makes the ratio test keep growing h
  // until it probes the parameter far outside the region where the local model holds
  // (e.g. an eta perturbed by several units), producing a degenerate/failed solve
  // that corrupts the shared solver state for every later finite difference.  Too
  // small: repeated shrinking drives h into the roundoff-dominated regime (central
  // optimum ~ eps^(1/3)).  Clamp both ends so the step stays sane.
  double l = 0, u = R_PosInf, rcur = NA_REAL;
  double hlast = h;

  arma::vec fp1(f0.size());
  arma::vec fm1(f0.size());
  gr.zeros(f0.n_elem); // avoid uninitialized read if no finite central step is found

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
      if (h >= hMax) break; // already at the cap and still growing: stop here
      h = nu*h;
      if (h > hMax) h = hMax; // clamp; probe once at hMax, then break if still growing
    } else if (l == 0) {
      if (h <= hMin) break; // already at the floor and still shrinking: stop here
      h = h/nu;
      if (h < hMin) h = hMin; // clamp; probe once at hMin, then break if still shrinking
    } else {
      h = (l + u)/2.0;
    }
  }
  return h;
}

// Column k of shi21Hessian() at the fixed step h; x[k] is put back exactly.
static void shi21HessColumn(shi21fn_type grad, arma::vec &x, arma::vec &gr0, int id,
                            int type, int k, double h, arma::vec &col) {
  double xk = x[k];
  x[k] += h;
  arma::vec grPH = grad(x, id);
  bool forwardFinite = grPH.is_finite();
  if (type == shi21HessForward && forwardFinite) {
    col = (grPH - gr0)/h;
    x[k] = xk;
    return;
  }
  x[k] -= 2*h;
  arma::vec grMH = grad(x, id);
  x[k] = xk;
  bool backwardFinite = grMH.is_finite();
  if (forwardFinite && backwardFinite) {
    // only reached for central: forward returned above
    col = (grPH - grMH)/(2.0*h);
  } else if (forwardFinite) {
    col = (grPH - gr0)/h;
  } else if (backwardFinite) {
    col = (gr0 - grMH)/h;
  }
}

arma::mat shi21Hessian(shi21fn_type grad, arma::vec &x, arma::vec &gr0, int id,
                       int type, double *hh, double ef, int maxiter,
                       double hMax, const double *hMin) {
  arma::mat H(x.n_elem, x.n_elem, arma::fill::zeros);
  if (type != shi21HessForward && type != shi21HessCentral) return H;
  arma::vec col(x.n_elem);
  for (int k = x.n_elem; k--;) {
    double h = hh[k];
    if (h <= 0) {
      double hMinK = (hMin == NULL) ? shi21hMinDefault : hMin[k];
      hh[k] = (type == shi21HessForward) ?
        shi21Forward(grad, x, h, gr0, col, id, k, ef, 1.5, 6.0, maxiter, hMax, hMinK) :
        shi21Central(grad, x, h, gr0, col, id, k, ef, 1.5, 4.5, 3.0, maxiter, hMax, hMinK);
      H.col(k) = col;
      continue;
    }
    col.zeros();
    shi21HessColumn(grad, x, gr0, id, type, k, h, col);
    H.col(k) = col;
  }
  return 0.5*(H + H.t());
}
