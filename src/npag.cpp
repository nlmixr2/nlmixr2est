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
};

// The NPAG adaptive-grid cycle (Yamada Alg 1).  Requires the FOCEi inner solve
// already set up (vaeInnerSetup_ / foceiSetup_), so npBuildPsiCore can fill Psi.
// Returns the final support points (eta space), their weights, the log-
// likelihood, cycle count, and a converged flag.
struct npagResult {
  arma::mat theta;   // nspp x neta support points (eta space)
  arma::vec lambda;  // nspp weights (sum 1)
  double objf;       // maximized log-likelihood
  int cycle;
  bool converged;
};

static npagResult npagRunCycle(const arma::vec& lower, const arma::vec& upper,
                               const npagCtl& ctl) {
  arma::mat theta = npSobolGrid(ctl.points, lower, upper);
  double eps = ctl.epsInit, f0 = -1e30, f1 = 0.0, lastObj = -1e30, objf = R_NegInf;
  int cycle = 0;
  bool converged = false;
  arma::vec lam;
  for (;;) {
    cycle++;
    if (cycle > 1) {
      theta = npExpandGrid(theta, eps, lower, upper, ctl.thetaD);
    }
    // estimation: Psi then Burke weights
    arma::mat psi;
    npBuildPsiCore(theta, ctl.cores, psi);
    double obj0 = 0.0;
    lam = npBurke(psi, &obj0);
    // condensation: weight threshold, then QR rank-revealing, then re-solve
    arma::uvec wk = npCondenseWeights(lam, ctl.ratio);
    theta = theta.rows(wk); psi = psi.cols(wk);
    arma::uvec qk = npCondenseQR(psi, ctl.qrTol);
    theta = theta.rows(qk); psi = psi.cols(qk);
    lam = npBurke(psi, &objf);
    // evaluation (Yamada Alg 1 convergence controller)
    if (std::fabs(lastObj - objf) <= ctl.thetaG && eps > ctl.thetaE) {
      eps /= 2.0;
      if (eps <= ctl.thetaE) {
        f1 = accu(log(psi * lam));
        if (std::fabs(f1 - f0) <= ctl.thetaF) { converged = true; break; }
        f0 = f1; eps = ctl.epsInit;
      }
    }
    if (cycle >= ctl.cycles) break;
    lastObj = objf;
  }
  npagResult r;
  r.theta = theta; r.lambda = lam; r.objf = objf; r.cycle = cycle;
  r.converged = converged;
  return r;
}

// M4: the driver is exercised in isolation through npagCycle_ (below).  Wiring it
// into the full nlmixr2 fit object (nlmixr2EnvSetup, tables, covariance) is a
// later milestone; until then est="npag" reports that.
void npagOuter(Environment e) {
  stop("NPAG estimation ('est=\"npag\"') is not yet wired to the fit object "
       "(the adaptive-grid cycle is available via npagCycle_)");
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
//' @return A list with \code{support} (support points, eta space; one per row),
//'   \code{weights}, \code{objf} (log-likelihood), \code{cycles}, and
//'   \code{converged}.
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List npagCycle_(arma::vec lower, arma::vec upper, int points = 2028,
                      int cycles = 100, int cores = 1) {
  npagCtl ctl;
  ctl.points = points; ctl.cycles = cycles; ctl.cores = cores;
  npagResult r = npagRunCycle(lower, upper, ctl);
  Rcpp::NumericMatrix support(r.theta.n_rows, r.theta.n_cols);
  std::copy(r.theta.begin(), r.theta.end(), support.begin());
  Rcpp::NumericVector weights(r.lambda.begin(), r.lambda.end());
  return Rcpp::List::create(
    Rcpp::Named("support") = support,
    Rcpp::Named("weights") = weights,
    Rcpp::Named("objf") = r.objf,
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
