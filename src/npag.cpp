// NPAG -- Nonparametric Adaptive Grid (Yamada 2021).  The estimation loop lives
// entirely here and reuses the FOCEI inner likelihood machinery in inner.cpp
// (via the np.h / imp.h numeric interfaces) to fill the Psi matrix of per-subject
// conditional likelihoods at each support point.  Called from foceiFitCpp_ in
// place of foceiOuter when est=="npag".
#include <RcppArmadillo.h>
#include "np.h"
#include "npCommon.h"

using namespace Rcpp;

// NPAG tuning constants (Yamada 2021): eps convergence floor / objf tolerance /
// F-candidate tolerance / minimum scaled support-point distance.
struct npagCtl {
  int points = 2028;    // initial Sobol grid size
  int cycles = 100;     // maximum cycles
  int cores = 1;
  double ratio = 1e-3;  // weight-threshold condensation (max*ratio)
  double qrTol = 1e-8;  // QR rank-revealing tolerance
  double epsInit = 0.2; // initial adaptive-grid fraction
  double thetaE = 1e-4; // eps floor
  double thetaG = 1e-4; // objf-change tolerance that halves eps
  double thetaF = 1e-2; // successive-F tolerance for the exit test
  double thetaD = 1e-4; // minimum scaled distance for daughter points
  bool gammaOptimize = false; // optimize the residual-error magnitude (gamma)
  double gammaInit = 1.0;     // initial gamma multiplier
  double gammaDelta = 0.1;    // initial gamma step fraction
};

// The NPAG adaptive-grid cycle (Yamada Alg 1).  Requires the FOCEi inner solve
// already set up (vaeInnerSetup_ / foceiSetup_), so npBuildPsiCore can fill Psi.
// Returns the final support points (eta space), their weights, the log-
// likelihood, cycle count, and a converged flag.
//
// objf is the nonparametric marginal log-likelihood
//   sum_i log(sum_k lambda_k p(y_i | support point k))
// using the inner conditional-density constant convention (and, in the gamma
// path, a per-subject log-sum-exp offset).  This is NOT the NONMEM/FOCEI OFV
// convention, so the resulting -2LL is not directly comparable to nlmixr2's
// FOCEI/SAEM/FOCE -2LL -- compare NPAG to NPAG (or to Pmetrics NPAG) instead.
struct npagResult {
  arma::mat theta;   // nspp x neta support points (eta space)
  arma::vec lambda;  // nspp weights (sum 1)
  double objf;       // maximized log-likelihood
  double gamma;      // final residual-error multiplier (1 if not optimized)
  int cycle;
  bool converged;
};

static npagResult npagRunCycle(const arma::vec& lower, const arma::vec& upper,
                               const npagCtl& ctl) {
  arma::mat theta = npSobolGrid(ctl.points, lower, upper);
  double eps = ctl.epsInit, f0 = -1e30, f1 = 0.0, lastObj = -1e30, objf = R_NegInf;
  double gamma = ctl.gammaInit, gDelta = ctl.gammaDelta;
  int cycle = 0;
  bool converged = false;
  arma::vec lam;
  arma::mat psi;
  for (;;) {
    cycle++;
    if (cycle > 1) {
      theta = npExpandGrid(theta, eps, lower, upper, ctl.thetaD);
    }
    // estimation: Psi at the current gamma, then Burke weights.  The gamma path
    // builds Psi through likInner0 (npBuildPsiCoreScaled) so censoring and
    // transform-both-sides are handled at the scaled residual error.
    double obj0 = 0.0, offset = 0.0;
    if (ctl.gammaOptimize) {
      npBuildPsiCoreScaled(theta, ctl.cores, gamma, psi, &offset);
    } else {
      npBuildPsiCore(theta, ctl.cores, psi);
    }
    lam = npBurke(psi, &obj0);
    // condensation: weight threshold, then QR rank-revealing
    arma::uvec wk = npCondenseWeights(lam, ctl.ratio);
    theta = theta.rows(wk); psi = psi.cols(wk);
    arma::uvec qk = npCondenseQR(psi, ctl.qrTol);
    theta = theta.rows(qk); psi = psi.cols(qk);
    // re-solve the weights on the condensed grid
    if (ctl.gammaOptimize) {
      npBuildPsiCoreScaled(theta, ctl.cores, gamma, psi, &offset);
      double bObj = 0.0; lam = npBurke(psi, &bObj); objf = bObj + offset;
    } else {
      lam = npBurke(psi, &objf);
    }
    // optimizations: residual-error magnitude (gamma) up/down search on the
    // condensed grid.  Objectives are offset-corrected (true log-likelihood).
    if (ctl.gammaOptimize) {
      double gUp = gamma * (1.0 + gDelta);
      double gDn = gamma / (1.0 + gDelta);
      arma::mat psiU, psiD;
      double offU = 0.0, offD = 0.0, bU = 0.0, bD = 0.0;
      npBuildPsiCoreScaled(theta, ctl.cores, gUp, psiU, &offU);
      arma::vec lamU = npBurke(psiU, &bU); double objU = bU + offU;
      npBuildPsiCoreScaled(theta, ctl.cores, gDn, psiD, &offD);
      arma::vec lamD = npBurke(psiD, &bD); double objD = bD + offD;
      if (objU > objf) { gamma = gUp; objf = objU; lam = lamU; psi = psiU; gDelta *= 4.0; }
      if (objD > objf) { gamma = gDn; objf = objD; lam = lamD; psi = psiD; gDelta *= 4.0; }
      gDelta *= 0.5;
      if (gDelta <= 0.01) gDelta = 0.1;
      if (gDelta > 1.0) gDelta = 1.0;         // cap the step so gamma cannot jump wildly
      gamma = std::min(std::max(gamma, 1e-3), 1e3);  // keep gamma in a sane range
    }
    // evaluation (Yamada Alg 1 convergence controller); objf is the true
    // (offset-corrected) log-likelihood in both paths.
    if (std::fabs(lastObj - objf) <= ctl.thetaG && eps > ctl.thetaE) {
      eps /= 2.0;
      if (eps <= ctl.thetaE) {
        f1 = objf;
        if (std::fabs(f1 - f0) <= ctl.thetaF) { converged = true; break; }
        f0 = f1; eps = ctl.epsInit;
      }
    }
    if (cycle >= ctl.cycles) break;
    lastObj = objf;
  }
  npagResult r;
  r.theta = theta; r.lambda = lam; r.objf = objf; r.gamma = gamma;
  r.cycle = cycle; r.converged = converged;
  return r;
}

// M4: the driver is exercised in isolation through npagCycle_ (below).  Wiring it
// into the full nlmixr2 fit object (nlmixr2EnvSetup, tables, covariance) is a
// later milestone; until then est="npag" reports that.
void npagOuter(Environment e) {
  stop("NPAG estimation ('est=\"npag\"') is not yet wired to the fit object "
       "(the adaptive-grid cycle is available via npagCycle_)");
}

//' Diagnostic: NPAG objective at a fixed grid and residual multiplier gamma
//' @param etaPoints support points, one per row
//' @param cores threads
//' @param gamma residual-error multiplier
//' @return offset-corrected marginal log-likelihood
//' @keywords internal
//' @export
// [[Rcpp::export]]
double npObjAtGamma_(arma::mat etaPoints, int cores, double gamma) {
  arma::mat psi; double offset = 0.0;
  npBuildPsiCoreScaled(etaPoints, cores, gamma, psi, &offset);
  double b = 0.0; arma::vec lam = npBurke(psi, &b);
  return b + offset;
}

//' Run the NPAG adaptive-grid cycle on a set-up inner problem
//'
//' Requires the FOCEi inner problem to be set up (\code{.npInnerSetup}).  Runs
//' the full Yamada adaptive-grid cycle (Sobol grid -> Psi -> Burke IPM ->
//' condensation -> expansion -> convergence) and returns the discrete mixing
//' distribution.  Exposed for testing ahead of the full fit-object wiring.
//'
//' @param lower,upper Numeric vectors, the per-eta support-point box.
//' @param points Initial Sobol grid size.
//' @param cycles Maximum cycles.
//' @param cores OpenMP threads.
//' @param gammaOptimize Optimize the residual-error magnitude (gamma) each cycle
//'   (only valid for uncensored normal endpoints).
//' @return A list with \code{support} (support points, eta space; one per row),
//'   \code{weights}, \code{objf} (log-likelihood), \code{gamma}, \code{cycles},
//'   and \code{converged}.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npagCycle_(arma::vec lower, arma::vec upper, int points = 2028,
                      int cycles = 100, int cores = 1, bool gammaOptimize = false) {
  npagCtl ctl;
  ctl.points = points; ctl.cycles = cycles; ctl.cores = cores;
  ctl.gammaOptimize = gammaOptimize;
  npagResult r = npagRunCycle(lower, upper, ctl);
  Rcpp::NumericMatrix support(r.theta.n_rows, r.theta.n_cols);
  std::copy(r.theta.begin(), r.theta.end(), support.begin());
  Rcpp::NumericVector weights(r.lambda.begin(), r.lambda.end());
  return Rcpp::List::create(
    Rcpp::Named("support") = support,
    Rcpp::Named("weights") = weights,
    Rcpp::Named("objf") = r.objf,
    Rcpp::Named("gamma") = r.gamma,
    Rcpp::Named("cycles") = r.cycle,
    Rcpp::Named("converged") = r.converged);
}

//' Burke interior-point weight solver (nonparametric maximum likelihood)
//'
//' Solves the convex nonparametric-maximum-likelihood weight problem for a fixed
//' set of support points: given the likelihood matrix \code{psi} (subjects in
//' rows, support points in columns) it returns the maximum-likelihood mixing
//' weights and the objective (log-likelihood).  Exposed for testing the C++
//' interior-point routine against golden fixtures.
//'
//' @param psi Numeric matrix, \code{psi[i, k] = p(y_i | support point k)}, with
//'   subjects in rows and support points in columns.
//' @return A list with \code{weights} (length \code{ncol(psi)}, non-negative,
//'   summing to 1) and \code{objective} (the maximized log-likelihood).
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npIpmBurke(arma::mat psi) {
  double obj = 0.0;
  arma::vec lam = npBurke(psi, &obj);
  Rcpp::NumericVector weights(lam.begin(), lam.end());  // plain vector, not Nx1
  return Rcpp::List::create(Rcpp::Named("weights") = weights,
                            Rcpp::Named("objective") = obj);
}

//' Sobol initial grid over a box (nonparametric engines)
//'
//' @param n Number of support points.
//' @param lower,upper Numeric vectors giving the per-dimension box bounds.
//' @return Numeric matrix, one support point per row.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix npSobolGrid_(int n, arma::vec lower, arma::vec upper) {
  arma::mat g = npSobolGrid(n, lower, upper);
  Rcpp::NumericMatrix out(g.n_rows, g.n_cols);
  std::copy(g.begin(), g.end(), out.begin());
  return out;
}

//' Condense support points (nonparametric engines)
//'
//' @param lambda Support-point weights.
//' @param psi Conditional-likelihood matrix (subjects x support points).
//' @param ratio Weight-threshold ratio (keep weight > max*ratio).
//' @param tol QR rank-revealing tolerance.
//' @return List with 1-based kept indices from the weight threshold
//'   (\code{weightKeep}) and from the subsequent QR pass (\code{qrKeep}).
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npCondense_(arma::vec lambda, arma::mat psi, double ratio = 1e-3,
                       double tol = 1e-8) {
  arma::uvec wk = npCondenseWeights(lambda, ratio);
  // apply the weight keep to psi columns, then QR-condense the survivors
  arma::mat psiW = psi.cols(wk);
  arma::uvec qk = npCondenseQR(psiW, tol);
  arma::uvec qkOrig = wk.elem(qk);   // map QR keep back to original indices
  return Rcpp::List::create(
    Rcpp::Named("weightKeep") = Rcpp::IntegerVector(wk.begin(), wk.end()) + 1,
    Rcpp::Named("qrKeep") = Rcpp::IntegerVector(qkOrig.begin(), qkOrig.end()) + 1);
}
