#include "nonMuThetaGrad.h"

void nonMuGradAccumObs(nonMuObjKind kind, double y, double f,
                       double gsd, double dgsdf,
                       const double *dfdth, int nth, double w,
                       arma::vec &score, arma::mat &info) {
  if (nth <= 0 || dfdth == NULL) return;
  if (!std::isfinite(f)) return;
  if (!std::isfinite(w) || w <= 0.0) return;
  for (int a = 0; a < nth; ++a) if (!std::isfinite(dfdth[a])) return;

  if (kind == nonMuObjLl) {
    // f IS the per-observation log-likelihood.  Minimizing -sum(f), so the
    // score is -d(f)/d(theta); the information is the BHHH outer product of the
    // per-observation score, which is what NONMEM's own non-mu route uses
    // (technical guide eqs. 1.47-1.52).
    for (int a = 0; a < nth; ++a) {
      score(a) -= w*dfdth[a];
      for (int b = 0; b < nth; ++b) info(a, b) += w*dfdth[a]*dfdth[b];
    }
    return;
  }

  if (!std::isfinite(y) || !std::isfinite(gsd) || !(gsd > 0.0)) return;
  const double r = y - f;
  const double g2 = gsd*gsd;
  // The residual scale moves with the prediction for a proportional or combined
  // error model, so d(g)/d(theta) = (dg/df)*(df/dtheta) and BOTH the mean and
  // the scale contribute.  Dropping the second term is what makes a naive
  // implementation quietly wrong on prop() models.
  const double dgf = std::isfinite(dgsdf) ? dgsdf : 0.0;
  // d/df of  0.5*(r/g)^2 + log(g)
  const double dObjdf = -r/g2 + (-(r*r)/(g2*gsd) + 1.0/gsd)*dgf;
  // Fisher information of a Gaussian in (mean, sd), both moving with theta:
  //   1/g^2 (df/dth)(df/dth)' + 2/g^2 (dg/dth)(dg/dth)'
  const double wMean = 1.0/g2;
  const double wScale = 2.0*dgf*dgf/g2;
  // extended least squares is twice the Gaussian -loglik; same derivatives
  const double sc = w*((kind == nonMuObjEls) ? 2.0 : 1.0);
  for (int a = 0; a < nth; ++a) {
    score(a) += sc*dObjdf*dfdth[a];
    for (int b = 0; b < nth; ++b) {
      info(a, b) += sc*(wMean + wScale)*dfdth[a]*dfdth[b];
    }
  }
}

bool nonMuGradStep(const arma::vec &score, const arma::mat &info,
                   const arma::vec &cur, double trust, arma::vec &step) {
  const arma::uword n = score.n_elem;
  if (n == 0 || info.n_rows != n || info.n_cols != n) return false;
  if (!score.is_finite() || !info.is_finite()) return false;
  if (cur.n_elem != n) return false;
  if (!(trust > 0.0)) trust = 0.5;

  // Relative to each theta's own magnitude, floored at 1, so the bound means
  // the same thing for a theta of 0.1 and one of 1e5.
  auto relMag = [&](const arma::vec &s) {
    double m = 0.0;
    for (arma::uword k = 0; k < n; ++k) {
      double d = std::fabs(cur(k));
      if (d < 1.0) d = 1.0;
      double q = std::fabs(s(k))/d;
      if (q > m) m = q;
    }
    return m;
  };

  double rc = 0.0;
  if (info.is_finite()) rc = arma::rcond(info);
  const bool wellCond = R_FINITE(rc) && rc > 1e-10;
  arma::vec cand;
  // Trustworthy information: take the exact step.  A large step off a
  // well-determined Hessian is a legitimate move, and damping it would slow
  // every healthy fit -- the sanity bound here is deliberately loose, meant
  // only to catch a step no honest refinement would propose.
  if (wellCond && arma::solve(cand, info, -score) && cand.is_finite() &&
      relMag(cand) <= 20.0) {
    step = cand;
    return true;
  }
  // Otherwise Levenberg-Marquardt, escalated until the step lands inside the
  // trust region.  The information is Gauss-Newton (or BHHH) and hence positive
  // semi-definite, so adding lambda*diag only improves conditioning and rotates
  // the step toward a scaled gradient step.
  double hscale = arma::abs(info.diag()).max();
  if (!R_FINITE(hscale) || hscale <= 0.0) hscale = 1.0;
  double lambda = 1e-8;
  for (int k = 0; k < 12; ++k) {
    arma::mat d = info;
    d.diag() += lambda*hscale;
    if (arma::solve(cand, d, -score) && cand.is_finite() &&
        relMag(cand) <= trust) {
      step = cand;
      return true;
    }
    lambda *= 100.0;
  }
  return false;
}
