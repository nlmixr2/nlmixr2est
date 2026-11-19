#ifndef __TRUNCNORM_H__
#define __TRUNCNORM_H__

// Standard-normal truncated draw, ported from rxode2's rxRmvn machinery
// (src/rxthreefry.cpp: ntail/tn/trandn), which is what backs
// truncnorm()/censResid.h's CWRES simulation.  That version takes a raw
// sitmo::threefry& engine, which nlmixr2est has no access to; here the SAME
// algorithm draws from rxNormEng()/rxUnifEng() instead -- the per-thread,
// threefry-seeded engine this file already uses everywhere else (see
// _saemSeedDoMcmc/_saemSeedCensAug), so this stays exactly as thread-safe
// and reproducible-under-seed as the rest of the SAEM chain, without a
// round-trip through R's SEXP/PROTECT machinery (unsafe from a non-main
// thread) that the rxRmvn R-facing entry point requires.
//
// Reference: Z. I. Botev (2015), "The Normal Law Under Linear Restrictions:
// Simulation and Estimation via Minimax Tilting", JRSS(B).  Chosen over a
// plain inverse-CDF draw (qnorm(pnorm(l) + u*(pnorm(u)-pnorm(l)))) because
// the inverse-CDF form loses precision catastrophically once l/u are a few
// SDs from 0 -- exactly the regime a BQL (M3/M4) observation's truncation
// region often sits in -- while trandn()'s accept-reject branches stay
// accurate arbitrarily far into the tail.

#include <cmath>

// Rayleigh-distribution accept-reject tail sampler (Marsaglia 1964),
// valid for a lower bound l > 0 (or, via trandn()'s sign flip, u < 0).
static inline double rxTruncNormNtail(double l, double u) {
  double c = l * l / 2.0;
  double f = expm1(c - u * u / 2.0);
  double x = 0.0;
  bool accept = true;
  while (accept) {
    double tmp = rxUnifEng(0.0, 1.0);
    x = c - log1p(rxUnifEng(0.0, 1.0) * f); // sample using Rayleigh
    if (tmp * tmp * x <= c) { // accepted
      accept = false;
    }
  }
  return sqrt(2.0 * x);
}

// Interior-region sampler: accept-reject from N(0,1) when the interval is
// wide enough that rejection is cheap; inverse-CDF (via R::qnorm/R::pnorm,
// safe here since the interval is narrow and near 0, so no cancellation)
// otherwise.
static inline double rxTruncNormTn(double l, double u, double tol = 2.05) {
  double x = 0.0;
  if (fabs(u - l) > tol) {
    x = rxNormEng(0.0, 1.0);
    while (x < l || x > u) {
      x = rxNormEng(0.0, 1.0);
    }
  } else {
    double pl = R::pnorm(l, 0.0, 1.0, true, false);
    double pu = R::pnorm(u, 0.0, 1.0, true, false);
    x = R::qnorm(pl + (pu - pl) * rxUnifEng(0.0, 1.0), 0.0, 1.0, true, false);
  }
  return x;
}

// Draw a standard normal truncated to [l, u] (+/-Inf accepted on either
// side).  a=0.4/tol=2.05 match rxode2's trandn() defaults exactly, so the
// method-selection behavior (and its accuracy characteristics) are identical
// to what CWRES's censored-observation simulation already relies on.
static inline double rxTruncNorm(double l, double u, double a = 0.4,
                                 double tol = 2.05) {
  if (l > a) {
    return rxTruncNormNtail(l, u);
  } else if (u < -a) {
    return -rxTruncNormNtail(-u, -l);
  }
  return rxTruncNormTn(l, u, tol);
}

#endif // __TRUNCNORM_H__
