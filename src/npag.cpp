// NPAG -- Nonparametric Adaptive Grid (Yamada 2021).  The estimation loop lives
// entirely here and reuses the FOCEI inner likelihood machinery in inner.cpp
// (via the np.h / imp.h numeric interfaces) to fill the Psi matrix of per-subject
// conditional likelihoods at each support point.  Called from foceiFitCpp_ in
// place of foceiOuter when est=="npag".
#include <RcppArmadillo.h>
#include "np.h"
#include "npCommon.h"

using namespace Rcpp;

// Module M0: scaffold only.  Later milestones assemble the adaptive-grid cycle
// (Psi -> Burke IPM -> condense -> gamma optimization -> expand -> converge).
void npagOuter(Environment e) {
  stop("NPAG estimation ('est=\"npag\"') is not yet implemented");
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
