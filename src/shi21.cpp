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
    gr = (f1-f0)/h;
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
      // Use last calculated gradient
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
  all = 1.0/all;
  double sum = 0.0;
  for (unsigned int j = all.size(); j--;) {
    sum += all[j];
  }
  return (((double)all.size())/sum);
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
    h = pow(3.0*ef, 0.3333333333333333333333); 
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

double shiRS(double &h, shi21fn_type f, double ef, arma::vec &t, int &id, int &idx,
             arma::vec &fp1, arma::vec &fm1,
             arma::vec &fp2, arma::vec &fm2, 
             double &l, double &u,
             bool &finiteFp1, bool &finiteFp2, bool &finiteFp4,
             bool &finiteFm1, bool &finiteFm2, bool &finiteFm4) {
  arma::vec tp4 = t;
  arma::vec tp2 = t;
  arma::vec tp1 = t;
  arma::vec tm4 = t;
  arma::vec tm2 = t;
  arma::vec tm1 = t;
  tp4(idx) += 4*h;
  tp2(idx) += 2*h;
  tp1(idx) += h;
  tp4(idx) -= 4*h;
  tp2(idx) -= 2*h;
  tp1(idx) -= h;

  fp1 = f(tp1, id);
  finiteFp1 = fp1.is_finite();
  if (!finiteFp1) {
    finiteFp2 = true;
    finiteFp4 = true;
    finiteFm1 = true;
    finiteFm2 = true;
    finiteFm4 = true;
    return -1.0;
  }
  fm1 = f(tm1, id);
  finiteFm1 = fm1.is_finite();
  if (!finiteFm1) {
    finiteFp2 = true;
    finiteFp4 = true;
    finiteFm2 = true;
    finiteFm4 = true;
    return -1.0;
  }
  fp2 = f(tp2, id);
  finiteFp2 = fp2.is_finite();
  if (!finiteFm1) {
    finiteFp4 = true;
    finiteFm2 = true;
    finiteFm4 = true;
    return -1.0;
  }
  fm2 = f(tm2, id);
  finiteFm2 = fm2.is_finite();
  if (!finiteFm2) {
    finiteFp4 = true;
    finiteFm4 = true;
    return -1.0;
  }
  arma::vec fp4 = f(tp4, id);
  finiteFp4 = fp2.is_finite();
  if (!finiteFp4) {
    finiteFm4 = true;
    return -1.0;
  }
  arma::vec fm4 = f(tm2, id);
  finiteFm4 = fp2.is_finite();
  if (!finiteFp4) {
    return -1.0;
  }
  // alpha = 2
  // s = (-2, -1, 1, 2)
  // w = (1/12, -2/3, 2/3, -1/12) = c(1/12, -8/12, 8/12, -1/12)
  // alpha*s = (-4, -2, 2 -4)
  // alpha*w = (4/12, -32/12, 32/12, -4/12)
  // derivative = (1*f(x-2h)-8*f(x-h)+8*f(x+h)-f(x+2h))/(12*h)
  // alpha*f(x+s*alpha) = 4*f(x-4h)-32*f(x-2h)+32*f(x+2h)-4*f(x+4h)
  // derivative - alpha*f(x+s*alpha)
  // -4*f(x-4h) + 33*f(x-2h)-8*f(x-h)+8*f(x+h)-33*f(x+2h)+4*f(x+4h)
  // ||w|| = 324
  // 324*12=3888
  arma::vec all = abs(-4*fm4 + 33*fm2 - 8*fm1 + 8*fp1 - 33*fp2 + 4*fp4)/(3888.0*ef);
  if (fm4.size() == 1) {
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

double shi21Stencil(shi21fn_type f, arma::vec &t, double &h,
                    arma::vec &f0, arma::vec &gr, int id, int idx,
                    double ef, double rl, double ru, double nu,
                    int maxiter) {
// exp((3/10)*log(5/2)+(1/10)*log(13))=
// 1.701282120455891888611
// 1.701282120455891888611*pow(er, 0.2)
  if (h == 0) {
    h = 1.701282120455891888611*pow(ef, 0.2); 
  }
  double h0=h, tmp = h;
  double l = 0, u = R_PosInf, rcur = NA_REAL;
  double hlast = h;

  arma::vec fp1(f0.size());
  arma::vec fm1(f0.size());
  
  arma::vec fp2(f0.size());
  arma::vec fm2(f0.size());

  int iter=0;
  bool finiteFp1 = true, finiteFp2 = true, finiteFp4 = true,
    finiteFm1=true, finiteFm2=true,  finiteFm4=true,
    calcGrad = false;
  while(true) {
    iter++;
    if (iter > maxiter) {
      h = hlast;
      break;
    }
    rcur = shiRS(h, f, ef, t, id, idx,
                 fp1, fm1, fp2, fm2, l, u,
                 finiteFp1, finiteFp2, finiteFp4,
                 finiteFm1, finiteFm2, finiteFm4);
    if (rcur == -1.0) {
      // pick new h and try to calculate a gradient if needed.
      if (!finiteFp1) {
        // hnew*4 = hold*0.5
        h = h*0.5/4.0;
        continue;
      } else if (!finiteFm1) {
        if (!calcGrad) {
          // forward difference
          calcGrad = true;
          gr = (fp1-f0)/h;
        }
        h = h*0.5/4.0;
        continue;
      } else if (!finiteFp2) {
        // hnew*4 = hold*1.5
        if (!calcGrad) {
          calcGrad = true;
          gr = (fp1-fm1)/(2*h);
        }
        h = h*1.5/4.0;
        continue;
      } else if (!finiteFm2) {
        if (!calcGrad) {
          calcGrad = true;
          gr = (-2.0*fm1 - 3.0*f0 + 6.0*fp1 - fp2)/(6.0*h);
        }
        h = h*1.5/4.0;
        continue;
      }
      // hnew*4 = hold*3.5
      h = h*3.5/4.0;
      if (!calcGrad) {
        calcGrad = true;
        gr = (fm2-8*fm1+8*fp1-fp2)/(12.0*h);
        hlast = h;
      }
      continue;
    } else {
      // derivative = (f(x-2h)-8*f(x-h)+8*f(x+h)-f(x+2h))/(12*h)
      // save good gradient and good h
      calcGrad = true;
      gr = (fm2-8*fm1+8*fp1-fp2)/(12.0*h);
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
